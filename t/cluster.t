use 5.008;                   # Require at least Perl version 5.10
use strict;                  # Must declare all variables before using them
use warnings;                # Emit helpful warnings
use Test::More;              # Testing module
use Test::LongString;        # Compare strings byte by byte
use Data::Section -setup;    # Have various DATA sections, allows for mock files
use lib 'lib';               # add 'lib' to @INC
use Bio::App::SELEX::RNAmotifAnalysis;
use File::Slurp qw( slurp );

use autodie;    # Automatically throw fatal exceptions for common unrecoverable
                #   errors (e.g. trying to open a non-existent file)

{   # Test 1
    my $fh            = fh_from('input');
    my $expected_aref = [
        [ 'abcdefghijklmnop',   3 ],
        [ 'erwoprhasdfasfd',    2 ],
        [ 'ponmlkjihgfdcba',    2 ],
        [ 'ponmlkjihgfedcba',   1 ],
        [ 'abcdeghijklmnop',    1 ],
        [ 'abcddefghijklmnopp', 1 ],
        [ 'ponmlkjihgfedcbb',   1 ],
    ];
    test_get_sequences_from( $fh, $expected_aref );
}

sub test_get_sequences_from {
    my $fh            = shift;
    my $expected_aref = shift;
    my $result_aref   = Bio::App::SELEX::RNAmotifAnalysis::get_sequences_from($fh, 'simple');
    is_deeply( $result_aref, $expected_aref, 'sequences correctly extracted' );
}

{   #Test 2
    my $expected_cluster    = 3;
    my $sample_cluster_href = sample_cluster();
    my $seq                 = [ 'ponmlkjihgfedcbd', 'X' ];
    my $result_cluster =
      Bio::App::SELEX::RNAmotifAnalysis::matching_cluster( 5, $sample_cluster_href, $seq,
        5 );
    is( $result_cluster, $expected_cluster,
        'found correct cluster for ' . $seq->[1] );
}

{   #Test 3
    my $expected_cluster_href = sample_cluster();
    my $input_fh              = fh_from('input');
    my $result_cluster_href   = Bio::App::SELEX::RNAmotifAnalysis::cluster(
        max_distance => 5,
        fh           => $input_fh,
        max_clusters => 5,
        file_type    => 'simple',
    );

    is_deeply( $result_cluster_href, $expected_cluster_href,
        'clusters correctly determined' );
}

{   #Test 4
    my $input_fh            = fh_from('input');
    my $result_cluster_href = Bio::App::SELEX::RNAmotifAnalysis::cluster(
        max_distance => 5,
        fh           => $input_fh,
        max_clusters => 2,
        file_type    => 'simple',
    );
    my $expected_cluster_href = sample_cluster1();

    is_deeply( $result_cluster_href, $expected_cluster_href,
        'limiting number of clusters works!' );
}

{   #Test 5
    my $input_fh            = fh_from('odd');
    my $seed_fh             = fh_from('seed');
    my $result_cluster_href = Bio::App::SELEX::RNAmotifAnalysis::cluster(
        max_distance => 5,
        fh           => $input_fh,
        seed_fh      => $seed_fh,
        max_clusters => 5,
        file_type    => 'simple',
    );
    my $expected_cluster_href = expected_odd_cluster();
    is_deeply( $result_cluster_href, $expected_cluster_href,
        'explicit seed clusters works!' );
}

{   # Test 6: Changing seed can change outcome (compare to Test 5 results)
    my $input_fh            = fh_from('seed');
    my $seed_fh             = fh_from('odd');
    my $result_cluster_href = Bio::App::SELEX::RNAmotifAnalysis::cluster(
        max_distance => 5,
        fh           => $input_fh,
        seed_fh      => $seed_fh,
        max_clusters => 5,
        file_type    => 'simple',
    );
    my $expected_cluster_href = expected_odd_as_seed_cluster();
    is_deeply( $result_cluster_href, $expected_cluster_href,
        'Different seed clusters works!' );
}

sub sample_cluster1 {
    return {
        1 => [
            [ 'abcdefghijklmnop',   3 ],
            [ 'abcdeghijklmnop',    1 ],
            [ 'abcddefghijklmnopp', 1 ],
        ],
        2 => [ [ 'erwoprhasdfasfd', 2 ], ],

    };
}

sub sample_cluster {
    return {
        %{ sample_cluster1() },
        3 => [
            [ 'ponmlkjihgfdcba',  2 ],
            [ 'ponmlkjihgfedcba', 1 ],
            [ 'ponmlkjihgfedcbb', 1 ],
        ],
    };
}

sub expected_odd_cluster {
    return {
        1 => [
            [ 'ABCDEFGHIJKLMNO', 1 ],
            [ 'BBCDEFGHIJKLMBB', 1 ],
            [ 'ABCDBBBBIJKLMNO', 1 ],
        ]
    };
}

sub expected_odd_as_seed_cluster {
    return {
        1 => [ [ 'BBCDEFGHIJKLMBB', 1 ], [ 'ABCDEFGHIJKLMNO', 1 ], ],
        2 => [ [ 'ABCDBBBBIJKLMNO', 1 ], ]
    };
}

{    # Test
    my $all_cluster_string;
    open( my $fh_all, '>', \$all_cluster_string );

    my $cluster1_string;
    open( my $fh1, '>', \$cluster1_string );

    my $cluster2_string;
    open( my $fh2, '>', \$cluster2_string );

    my $cluster3_string;
    open( my $fh3, '>', \$cluster3_string );

    my $cluster_href = sample_cluster();
    my $fh_href      = {
        1 => $fh1,
        2 => $fh2,
        3 => $fh3,
    };

    Bio::App::SELEX::RNAmotifAnalysis::write_out_clusters(
        cluster_href => $cluster_href,
        fh_all       => $fh_all,
        fh_href      => $fh_href,
        max_top_seqs => 1000
    );
    my $expected_cluster1 = string_from('cluster1');
    my $expected_cluster2 = string_from('cluster2');
    my $expected_cluster3 = string_from('cluster3');
    is_string( $cluster1_string, $expected_cluster1, 'cluster1 good' );
    is_string( $cluster2_string, $expected_cluster2, 'cluster2 good' );
    is_string( $cluster3_string, $expected_cluster3, 'cluster3 good' );
}

{    # Test 7: Maximum top sequences
    my $all_cluster_string;
    open( my $fh_all, '>', \$all_cluster_string );

    my %string_for = (
        1 => '',
        2 => '',
        3 => '',
    );
    open( my $fh1, '>', \$string_for{1} );

    open( my $fh2, '>', \$string_for{2} );

    open( my $fh3, '>', \$string_for{3} );

    my $cluster_href = sample_cluster();
    my $fh_href      = {
        1 => $fh1,
        2 => $fh2,
        3 => $fh3,
    };

    my $input_fh = fh_from('input');

    Bio::App::SELEX::RNAmotifAnalysis::write_out_clusters(
        cluster_href => $cluster_href,
        fh_all       => $fh_all,
        fh_href      => $fh_href,
        max_top_seqs => 2
    );

    for my $cluster ( 1 .. 2 ) {
        my $expected_fasta_top = string_from("fasta_top_cluster_$cluster");
        my $expected_fasta_overage =
          string_from("fasta_overage_cluster_$cluster");
        my $overage_file_name = "cluster_${cluster}_overage.fasta";
        my $result_overage    = slurp $overage_file_name;

        is( $string_for{$cluster}, $expected_fasta_top, 'top seqs works' );
        is( $result_overage, $expected_fasta_overage,
            'overage fasta file correct' );
        delete_temp_file($overage_file_name);
    }
}

{    # Test
    my $all_cluster_string;
    open( my $fh_all, '>', \$all_cluster_string );

    my $cluster1_string;
    open( my $fh1, '>', \$cluster1_string );

    my $cluster2_string;
    open( my $fh2, '>', \$cluster2_string );

    my $cluster3_string;
    open( my $fh3, '>', \$cluster3_string );

    my $cluster_href = sample_cluster();
    my $fh_href      = {
        1 => $fh1,
        2 => $fh2,
        3 => $fh3,
    };

    Bio::App::SELEX::RNAmotifAnalysis::write_out_clusters(
        cluster_href => $cluster_href,
        fh_all       => $fh_all,
        fh_href      => $fh_href,
        max_top_seqs => 1000
    );
    my $expected_cluster1 = string_from('cluster1');
    my $expected_cluster2 = string_from('cluster2');
    my $expected_cluster3 = string_from('cluster3');
    is_string( $cluster1_string, $expected_cluster1, 'cluster1 good' );
    is_string( $cluster2_string, $expected_cluster2, 'cluster2 good' );
    is_string( $cluster3_string, $expected_cluster3, 'cluster3 good' );
}

done_testing();

sub delete_temp_file {
    my $filename = shift;
    my $result   = unlink $filename;
    ok( $result, "successfully deleted temporary file '$filename'" );
}

sub fh_to_empty_string {
    my $string = '';
    open( my $fh, '>', \$string );
    return $fh;
}

sub sref_from {
    my $section = shift;

    #Scalar reference from the section
    return __PACKAGE__->section_data($section);
}

sub string_from {
    my $section = shift;

    #Get the scalar reference
    my $sref = sref_from($section);

    #Return the actual scalar (probably a string), not the reference to it
    return ${$sref};
}

sub fh_from {
    my $section = shift;
    my $sref    = sref_from($section);

    #Create filehandle to the referenced scalar
    open( my $fh, '<', $sref );
    return $fh;
}

#------------------------------------------------------------------------
# IMPORTANT!
#
# Each line from each section automatically ends with a newline character
#------------------------------------------------------------------------

__DATA__
__[ input ]__
abcdefghijklmnop
abcdefghijklmnop
abcdefghijklmnop
abcdeghijklmnop
abcddefghijklmnopp
ponmlkjihgfedcba
ponmlkjihgfedcbb
ponmlkjihgfdcba
ponmlkjihgfdcba
erwoprhasdfasfd
erwoprhasdfasfd
__[ cluster1 ]__
>1.1.3
abcdefghijklmnop
>1.2.1
abcdeghijklmnop
>1.3.1
abcddefghijklmnopp
__[ cluster2 ]__
>2.1.2
ponmlkjihgfdcba
>2.2.1
ponmlkjihgfedcba
>2.3.1
ponmlkjihgfedcbb
__[ cluster3 ]__
>3.1.2
erwoprhasdfasfd
>3.1.2b
erwoprhasdfasfd
__[ seed ]__
ABCDEFGHIJKLMNO
__[ odd ]__
BBCDEFGHIJKLMBB
ABCDBBBBIJKLMNO
__[ odd_cluster ]__
>1.1.1
ABCDBBBBIJKLMNO
>1.2.1
ABCDEFGHIJKLMNO
>1.3.1
BBCDEFGHIJKLMBB
__[ seeded_cluster1 ]__
>1.1.1
abcdeghijklmnop
>1.2.3
abcdefghijklmnop
>1.3.1
abcddefghijklmnopp
__[ fasta_top_cluster_1 ]__
>1.1.3
abcdefghijklmnop
>1.2.1
abcdeghijklmnop
__[ fasta_overage_cluster_1 ]__
>1.3.1
abcddefghijklmnopp
__[ fasta_top_cluster_2 ]__
>2.1.2
ponmlkjihgfdcba
>2.2.1
ponmlkjihgfedcba
__[ fasta_overage_cluster_2 ]__
>2.3.1
ponmlkjihgfedcbb
