#!/usr/bin/env perl

use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib", "$FindBin::Bin/PerlLib");
use DBI;
use Sqlite_connect;
use CanvasXpress::Scatter2D;
use CanvasXpress::PlotOnLoader;
use URI::Escape;

use Data::Dumper;

$|++;

our $SEE = 0;

main: {
    
    my $cgi = new CGI();
    print $cgi->header();
    
    my %params = $cgi->Vars();

    if ($params{SEE}) {
        $SEE = 1;
    }
    
    my $sqlite_db = $params{sqlite} or die "Error, need sqlite param";
        
    my $sample_pair = $params{pair} or die "Error, need sample pair";
    #$sample_pair =~ s/\.(genes|isoforms)\.results//g;
    
    my ($sampleA, $sampleB) = split(/,/, $sample_pair);

    my $plot_loader_func_name = "load_plots_$$";

    my $plot_loader = new CanvasXpress::PlotOnLoader($plot_loader_func_name);
    

    print $cgi->start_html(-title => "Comparison of sample $sampleA to $sampleB",
                           -onLoad => $plot_loader_func_name . "();",
                           -style => {'src' => "CSS/common.css"},
        );    
    
    
    my $feature_type = $params{feature_type} || 'G'; # default to Gene
    
    my ($genes_selected, $trans_selected, $both_selected) = ("","","");
    
    if ($feature_type eq 'G') {
        $genes_selected = "checked";
    }
    elsif ($feature_type eq 'T') {
        $trans_selected = "checked";
    }
    elsif ($feature_type eq 'B') {
        $both_selected = "checked";
        $feature_type = undef; # required to pull both features in API
    }

    
    my $dbproc = DBI->connect( "dbi:SQLite:$sqlite_db" ) || die "Cannot connect: $DBI::errstr";

    
    my %gene_to_data = &get_diff_express_data($dbproc, $sampleA, $sampleB, $feature_type);
    
    # print "<pre>" . Dumper(\%gene_to_data) . "</pre>\n";

    my $form_text = "<form action='DE_sample_pair.cgi' method='get'>\n"
                  . "  Show:  <input type='radio' name='feature_type' value='G' $genes_selected>Genes</input>\n"
                  . "         <input type='radio' name='feature_type' value='T' $trans_selected>Transcripts</input>\n"
                  . "         <input type='radio' name='feature_type' value='B' $both_selected>Both</input>\n"
                  . "     <input type='submit' value='Go'/>\n"
                  . "     <input type='hidden' name='sqlite' value=\'$sqlite_db\'/>\n"
                  . "     <input type='hidden' name='pair' value=\'$sample_pair\'/>\n"
                  . "</form>\n";
    
    print $form_text;
        
    print "<h2>MA plot: $sampleA vs. $sampleB</h2>\n";
    my $ma_plot_obj = &write_MA_plot($sampleA, $sampleB, \%gene_to_data, $sqlite_db);
    $plot_loader->add_plot($ma_plot_obj);

    print "<div style='clear:both;' />\n";
    print "<h2>Volcano plot: $sampleA vs. $sampleB</h2>\n";
    my $volcano_plot = &write_Volcano_plot($sampleA, $sampleB, \%gene_to_data, $sqlite_db);
    $plot_loader->add_plot($volcano_plot);
    
    print $plot_loader->write_plot_loader();
    
    print $cgi->end_html();
    
    exit(0);
}


####
sub get_diff_express_data {
    my ($dbproc, $sampleA, $sampleB, $feature_type) = @_;
    
    my %comparisons = &get_sample_comparisons_list($dbproc);  # sampleA . $; . sampleB
    
    my $invert_FC = 0;
    my ($sample_id_A, $sample_id_B); # just for querying the database below
    if ($comparisons{ join("$;", $sampleA, $sampleB) }) {
        ($sample_id_A, $sample_id_B) = ($sampleA, $sampleB);
    }
    elsif ($comparisons{ join("$;", $sampleB, $sampleA) }) {
        $invert_FC = 1;
        ($sample_id_A, $sample_id_B) = ($sampleB, $sampleA);
    }
    else {
        die "Error, did not detect an analysis of DE between $sampleA and $sampleB. Only have record of: " . Dumper(\%comparisons);
    }
    

   
    my $query = "select d.feature_name, d.log_avg_expr, d.log_fold_change, d.p_value, d.fdr, T.annotation "
        . " from Diff_expression d, Samples s1, Samples s2, Transcript T "
        . " where d.sample_id_A = s1.sample_id "
        . " and d.sample_id_B = s2.sample_id "
        . " and s1.sample_name = \"$sample_id_A\" "
        . " and s2.sample_name = \"$sample_id_B\" ";
    
    if ($feature_type) {
        
        $query .= " and d.feature_type = \"$feature_type\" ";
        $query .= ($feature_type eq 'G') ? " and T.gene_id = d.feature_name " : " and T.transcript_id = d.feature_name";
        
    }
    else {
        $query .= " and d.feature_name in (T.gene_id, T.transcript_id) ";
    }
    
    print STDERR "$query\n";
    my $start = time();
    # exit(0);

    my @results = &do_sql_2D($dbproc, $query);
    my $end = time();
    
    my $query_time = $end - $start;
    print "<p>Query time: $query_time seconds.</p>";
        

    my %gene_info;
    my %gene_annot;

    foreach my $result (@results) {
        my ($feature_name, $log_avg_expr, $log_fold_change, $p_value, $fdr, $annot) = @$result;

        if ($invert_FC) {
            $log_fold_change *= -1; # logspace
        }
        
        $gene_info{$feature_name} = { log_avg_expr => $log_avg_expr,
                                      log_fold_change => $log_fold_change,
                                      p_value => $p_value,
                                      fdr => $fdr,
                                      annot => $annot,
                                      
        };

        
    }
    
    return(%gene_info);
                                          

}



####
sub get_sample_comparisons_list {
    my ($dbproc) = @_;
        
    ## get relevant data:
    my $query = "select distinct s1.sample_name, s2.sample_name "
        . " from Samples s1, Samples s2, Diff_expression d "
        . " where s1.sample_id = d.sample_id_A "
        . " and s2.sample_id = d.sample_id_B ";
    my @results = &do_sql_2D($dbproc, $query);
    
    my %pairs;
    foreach my $result (@results) {
        my ($sampleA, $sampleB) = @$result;
        # print "<p>$sampleA, $sampleB\n";
        my $token = join("$;", $sampleA, $sampleB);
        $pairs{$token} = 1;
    }

    return(%pairs);

}


####
sub write_MA_plot {
    my ($sampleA, $sampleB, $gene_to_data_href, $sqlite_db) = @_;
    

    my @value_matrix;
    my @annots;
    my @stat_signif;
    foreach my $feature_name (keys %$gene_to_data_href) {
        my $struct = $gene_to_data_href->{$feature_name};
        
        my $log_avg_expr = $struct->{log_avg_expr};
        my $log_fold_change = $struct->{log_fold_change};
        my $fdr = $struct->{fdr};
        
        $log_avg_expr = sprintf("%.2f", $log_avg_expr);
        $log_fold_change = sprintf("%.2f", $log_fold_change);
        

        push (@value_matrix, [$feature_name, $log_avg_expr, $log_fold_change]);
        
        my $annot = $struct->{annot};
        push (@annots, $annot);
        
        my $is_stat_signif = ($fdr <= 0.05) ? "Yes" : "No";
        push (@stat_signif, $is_stat_signif);
    }

    $sqlite_db = uri_escape($sqlite_db);
    
    
    my %plot_inputs = ( replicate_names => ['log_avg_expr', 'log_fold_change'],
                        value_matrix => \@value_matrix,
                        comparisons => [ ['log_avg_expr', 'log_fold_change'] ],
                        
                        feature_annotations => {
                            annotation => \@annots,
                            significant => \@stat_signif,
                        },

                        colorBy => 'significant',
                        
                            #events => {

                            #'click' => "console.log('clicked');\n"
                            #    . "var gene = o['y']['vars'][0];\n"
                            #    . "console.log(gene);\n"
                            #    . "alert('clicked');\n"
                            #    ,
                            #    
                            #    'dblclick' => "alert('dblclick');\n",
                            #
                            #    'mousemove' => #"DumperAlert(o);",
                                #"var gene = o['vars'][0];\n"
                            #    "console.log('1');\n"
                            #    . "alert('mouseover');\n",
                        
                        

                        events =>  { 'dblclick' => "var gene = o['y']['vars'][0];\n"
                                         #. "document.location.href=\'feature_report.cgi?feature_name=\' + gene"
                                         #. " + \'&sqlite=$sqlite_db\';\n",
                                         . "window.open(\'feature_report.cgi?feature_name=\' + gene + \'&sqlite=$sqlite_db\');\n", 
                                         
                        },
                        
                        
                        
        );
    
    my $plot_obj = new CanvasXpress::Scatter2D("ma_plot_$$");

    print $plot_obj->draw(%plot_inputs);
    
    
    return ($plot_obj);
}


####
sub write_Volcano_plot {
    my ($sampleA, $sampleB, $gene_to_data_href, $sqlite_db) = @_;
    
    my @value_matrix;
    my @annots;
    my @stat_signif;
    foreach my $feature_name (keys %$gene_to_data_href) {
        my $struct = $gene_to_data_href->{$feature_name};
        
        my $fdr = $struct->{fdr};
        my $log_fold_change = $struct->{log_fold_change};
        
        my $is_stat_signif = ($fdr <= 0.05) ? "Yes" : "No";
        push (@stat_signif, $is_stat_signif);
        
        $fdr = -1 * log($fdr)/log(10);

        $log_fold_change = sprintf("%.2f", $log_fold_change);
        $fdr = sprintf("%.2f", $fdr);
        
        push (@value_matrix, [$feature_name, $log_fold_change, $fdr]);
        
        my $annot = $struct->{annot};
        push (@annots, $annot);
        

    }

    $sqlite_db = uri_escape($sqlite_db);
    
    my %plot_inputs = ( replicate_names => ['log_fold_change', '-1*log10(fdr)'],
                        value_matrix => \@value_matrix,
                        comparisons => [ ['log_fold_change', '-1*log10(fdr)'] ],
    
                        feature_annotations => {
                            annotation => \@annots,
                            significant => \@stat_signif,
                        },

                        colorBy => 'significant',


                        
                        events =>  { 'dblclick' => "var gene = o['y']['vars'][0];\n"
                                         #. "document.location.href=\'feature_report.cgi?feature_name=\' + gene"
                                         #. " + \'&sqlite=$sqlite_db\';\n",
                                         . "window.open(\'feature_report.cgi?feature_name=\' + gene + \'&sqlite=$sqlite_db\');\n", 

                        },

        );
    
    my $plot_obj = new CanvasXpress::Scatter2D("volcano_plot_$$");

    print $plot_obj->draw(%plot_inputs);
    
    
    return ($plot_obj);
}
