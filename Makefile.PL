
use strict;
use warnings;



use ExtUtils::MakeMaker 6.30;



my %WriteMakefileArgs = (
  "ABSTRACT" => "Cluster SELEX sequences and calculate their structures",
  "AUTHOR" => "Christopher Bottoms <molecules\@cpan.org>",
  "BUILD_REQUIRES" => {
    "Data::Section" => 0,
    "File::Slurp" => 0,
    "Test::LongString" => 0,
    "Test::More" => "0.96"
  },
  "CONFIGURE_REQUIRES" => {
    "ExtUtils::MakeMaker" => "6.30"
  },
  "DISTNAME" => "Bio-App-SELEX-RNAmotifAnalysis",
  "EXE_FILES" => [
    "bin/selex_covarianceSearch",
    "bin/selex_adapter_trimmer",
    "bin/RNAmotifAnalysis",
    "bin/selex_clustal2stockholm.pl",
    "bin/selex_stock2fasta.pl"
  ],
  "LICENSE" => "perl",
  "NAME" => "Bio::App::SELEX::RNAmotifAnalysis",
  "PREREQ_PM" => {
    "Config::Tiny" => 0,
    "Getopt::Long::Descriptive" => 0,
    "List::MoreUtils" => 0,
    "Text::LevenshteinXS" => 0,
    "autodie" => "2.0"
  },
  "VERSION" => "0.0084",
  "test" => {
    "TESTS" => "t/*.t"
  }
);


unless ( eval { ExtUtils::MakeMaker->VERSION(6.56) } ) {
  my $br = delete $WriteMakefileArgs{BUILD_REQUIRES};
  my $pp = $WriteMakefileArgs{PREREQ_PM};
  for my $mod ( keys %$br ) {
    if ( exists $pp->{$mod} ) {
      $pp->{$mod} = $br->{$mod} if $br->{$mod} > $pp->{$mod};
    }
    else {
      $pp->{$mod} = $br->{$mod};
    }
  }
}

delete $WriteMakefileArgs{CONFIGURE_REQUIRES}
  unless eval { ExtUtils::MakeMaker->VERSION(6.52) };

WriteMakefile(%WriteMakefileArgs);



