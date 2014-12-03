#!/usr/bin/env perl

use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib", "$FindBin::Bin/PerlLib");
use DBI;
use Sqlite_connect;
use Trinotate;
use URI::Escape;
use Data::Dumper;


our $SEE = 0;

my $MAX_RESULTS = 200;
my $MAX_TEXTLINE_LEN = 300;

main: {
    
    my $cgi = new CGI();
    print $cgi->header();
    
    my %params = $cgi->Vars();
    
    my $sqlite_db = $params{sqlite} or die "Error, need sqlite param";
    my $keyword = $params{keyword} or die "Error, need keyword param";
    $keyword =~ s/\W/ /g;
    
    print $cgi->start_html(-title => 'Keyword Search: $keyword',
                           -style => {'src' => "CSS/common.css"},
        );
    
    my $dbproc = DBI->connect( "dbi:SQLite:$sqlite_db" ) || die "Cannot connect: $DBI::errstr";

    my $query = "select gene_id, transcript_id, annotation from Transcript where ";
    my @keywords = split(/\s+/, $keyword);
    for (my $i = 0; $i <= $#keywords; $i++) {
        $query .= " annotation like \"%" . $keywords[$i] . "%\"";
        if ($i != $#keywords) {
            $query .= " and ";
        }
    }


    print "<h1>Search results for [$keyword]</h1>\n";
    
    my $counter = 0;

    my @results = &do_sql_2D($dbproc, $query);

    print "<p>There are " . scalar(@results) . " matching entries.</p>\n";
    
    if (@results) {
        
        print "<table border=1 >\n";
        print "<tr><th>#</th><th>gene_id</th><th>transcript_id</th><th>annotation</th></tr>\n";
        
        
        
        foreach my $result (@results) {
            
            $counter++;

            my ($gene_id, $transcript_id, $annot) = @$result;
            if (length($annot) > $MAX_TEXTLINE_LEN) {
                $annot = substr($annot, 0, $MAX_TEXTLINE_LEN);
            }
            foreach my $keyw (@keywords) {
                $annot =~ s/$keyw/<b>$keyw<\/b>/g;
            }
            print "<tr><td>$counter</td>"
                . "<td><a href=\"feature_report.cgi?feature_name=" . uri_escape($gene_id) . "&sqlite=" . uri_escape($sqlite_db) . "\" target=_blank >$gene_id</td>"
                . "<td><a href=\"feature_report.cgi?feature_name=" . uri_escape($transcript_id) . "&sqlite=" . uri_escape($sqlite_db) . "\" target=_blank >$transcript_id</td>"
                . "<td>$annot</td></tr>\n";
                    
            if ($counter >= $MAX_RESULTS) {
                print "<tr><td colspan=4>RESULTS TRUNCATED TO MAX OF $MAX_RESULTS ENTRIES</td></tr>\n";
                last;
            }
        }
        print "</table>\n";
    }
    
        
    $dbproc->disconnect;
        
    print $cgi->end_html();
    
    exit(0);

}

