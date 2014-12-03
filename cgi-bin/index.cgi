#!/usr/bin/env perl

use strict;
use warnings;

# standard perl modules
use CGI;
use CGI::Pretty ":standard";
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use DBI;
use URI::Escape;
use Data::Dumper;
use Cwd;

# custom modules
use lib ("$FindBin::Bin/../../PerlLib", "$FindBin::Bin/PerlLib");
use Sqlite_connect;


main: {
    
    my $cgi = new CGI();
    print $cgi->header();
    
    my %params = $cgi->Vars();

    print $cgi->start_html(-title => "Trinotate Web for Annotation and Expression Analysis (ultra-early release)",
                           -style => {'src' => "CSS/common.css"},
        );
    
    print "<h1>Trinotate Web for Annotation and Expression Analysis</h1>\n";
    

    my $sqlite_db = $params{sqlite_db};
    unless ($sqlite_db) {
        my $db = "TrinityFunctional.db";
        if (-s $db) {
            $sqlite_db = "$db";
        }
        
    }
                
    if ($sqlite_db) {
        
        if (! -s $sqlite_db) {
            die "Error, cannot locate $sqlite_db";
        }
        &TrinotateWebMain(\%params, $sqlite_db);
        
    }
    else {
        
        print "<h2>Need database info</h2>\n";
        
        &print_get_sqlite_db_path_form();
        
    }
        
    print $cgi->end_html();

    exit(0);
}

####
sub print_get_sqlite_db_path_form {
    
    print "<form action='index.cgi' method='get'>\n";

    print "<ul>\n"
        . "<li>Path to Trinotate SQLite database:\n"
        . "<li><input type='text' name='sqlite_db' maxlength=1000 size=150 />\n"
        . "<li><input type='submit' />\n"
        . "</li>\n"
        . "</ul>\n";
    
    print "</form>\n";


    return;
}
    
####
sub TrinotateWebMain {
    my ($params_href, $sqlite_db) = @_;
    
    my $dbproc = DBI->connect( "dbi:SQLite:$sqlite_db" ) || die "Cannot connect: $DBI::errstr";
        
    
    print "<h2>Stats</h2>\n" . &stats_panel_text($dbproc);
    
    print "<h2>Annotation Keyword Search</h2>\n" . &keyword_search_panel_text($dbproc, $sqlite_db);
    
    print "<h2>Gene or Transcript ID Search</h2>\n" . &gene_search_panel($dbproc, $sqlite_db);
    
    print "<h2>Pairwise Expression Comparisons (Volcano and MA plots)</h2>\n" . &pairwise_DE_panel_text($dbproc, $sqlite_db);

    print "<h2>Multi-sample Comparisons (Expression Profiling)</h2>\n" . &multi_sample_DE_panel_text($dbproc, $sqlite_db);

    
    return;
}

####
sub stats_panel_text {
    my ($dbproc) = @_;
    
    #my $query = "select count(distinct gene_id) from Transcript";
    #my $gene_count = &very_first_result_sql($dbproc, $query);
    
    #$query = "select count(distinct transcript_id) from Transcript";
    #my $transcript_count = &very_first_result_sql($dbproc, $query);

    my $html = "<p>Various summary stats go here...</p>\n";
    
    #$html .= "<p>Got $gene_count genes and $transcript_count transcripts</p>\n";
    
    return($html);
    
}

####
sub keyword_search_panel_text {
    my ($dbproc, $sqlite_db) = @_;
    
    my $html = "<p>Text search of transcript annotations</p>\n";

    $html .= "<form action='keyword_search.cgi'>\n"
        . "<input type='text' size=50 name='keyword' />\n"
        . "<input type='hidden' name='sqlite' value=\'$sqlite_db\' />\n"
        . "<input type='submit' value='search' />\n"
        . "</form>\n";
    
    $html .= "<p>Still needed: search based on specific attribute: pfam, go, kegg, etc.</p>";

    return($html);
}

####
sub gene_search_panel {
    my ($dbproc, $sqlite_db) = @_;
    
    my $html = "<p>Gene or Transcript Identifier to Feature Report:</p>\n";

    $html .= "<form action='feature_report.cgi'>\n"
        . "<input type='text' size=50 name='feature_name' />\n"
        . "<input type='hidden' name='sqlite' value=\'$sqlite_db\' />\n"
        . "<input type='submit' value='search' />\n"
        . "</form>\n";


    return($html);
}


####
sub pairwise_DE_panel_text {
    my ($dbproc, $sqlite_db) = @_;
    
    ## get list of pairs:
    my $query = "select sample_name from Samples";
    my @samples = &do_sql($dbproc, $query);
    ## generate the matrix
    my $html = "<table border=1 >\n";
    $html .= "<tr><th></th>";
    for (my $i = 1; $i <= $#samples; $i++) {
        $html .= "<th>$samples[$i]</th>";
    }
    $html .= "</tr>\n";
    
    for (my $i = 0; $i < $#samples; $i++) {
        $html .= "<tr><th>$samples[$i]</th>";
        
        my $sample_i = $samples[$i];

        for (my $j = 1; $j <= $#samples; $j++) {
            if ($j > $i) {
                my $sample_j = $samples[$j];
                my $compare_pair = "$sample_i vs. $sample_j";
                $html .= "<td>"
                    
                    ## link to volcano and MA plots for sample pair
                    . "<a href=\"DE_sample_pair.cgi?pair=" . uri_escape("$sample_i,$sample_j") . "&sqlite=" . uri_escape("$sqlite_db") . "\" "
                    . " target='_blank' >"
                    #. "$compare_pair"
                    . "[MA,Vo]"
                    . "</a>"

                    ## link to heatmap for the DE features for sample pair:
                    . "<a href=\"HeatmapNav.cgi?sqlite=" . uri_escape("$sqlite_db") . "&sample_pair=" . uri_escape("$sample_i,$sample_j") . "\" "
                    . " target='_blank' >"
                    . " [HMap]"
                    . " </a>"
                    . "</td>";
            }
            else {
                $html .= "<td></td>\n";
            }
        }
        $html .= "</tr>\n";
    }
    $html .= "</table>\n";
    
    return($html);
}


####
sub multi_sample_DE_panel_text {
    my ($dbproc, $sqlite_db) = @_;
    
    my $url = "HeatmapNav.cgi?sqlite=" . uri_escape("$sqlite_db");
    my $html = "<p>Go to the interactive <a href=\"$url\" target='_blank'>heatmap</a> for all DE transcripts.</p>\n"; #?sqlite=$sqlite_db\">heatmap navigator</a></p>\n";
    
    ## get cluster analyses
    my $query = "select ECA.cluster_analysis_id, ECA.cluster_analysis_group, ECA.cluster_analysis_name, count(distinct EC.expr_cluster_id) "
        . " from ExprClusterAnalyses ECA, ExprClusters EC "
        . " where ECA.cluster_analysis_id = EC.cluster_analysis_id "
        . " group by ECA.cluster_analysis_id, ECA.cluster_analysis_group, ECA.cluster_analysis_name order by ECA.cluster_analysis_group";
    
    #print $query;

    my @results = &do_sql_2D($dbproc, $query);
    if (@results) {
        $html .= "<ul>Analyses of clusters of expression profiles:\n";
        
        foreach my $result (@results) {
            my ($cluster_analysis_id, $analysis_group, $analysis_name, $cluster_count) = @$result;
            $html .= "<li><a href=\"transcript_cluster_viewer.cgi?cluster_analysis_id=" . uri_escape($cluster_analysis_id)
                . "&sqlite=" . uri_escape($sqlite_db) . "\" target=_blank >$analysis_group :: $analysis_name</a> with $cluster_count clusters.\n";
        }
        $html .= "</ul>\n";
    }
    else {
        $html .= "<p>No expression profile clusters defined yet.\n";
    }
    
    
    return($html);
}
