#!/usr/bin/env perl

use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib", "$FindBin::Bin/PerlLib");
use DBI;
use Sqlite_connect;
use CanvasXpress::Heatmap;
use CanvasXpress::PlotOnLoader;
use Data::Dumper;
use File::Basename;
use DEWebCommon;
use Data::Dumper;

$|++;

main: {
    
    my $cgi = new CGI();
    print $cgi->header();
    
    my %params = $cgi->Vars();
    my $sqlite_db = $params{sqlite} or die "Error, need sqlite param";

    my $sample_pair = $params{sample_pair};

    my $title = ($sample_pair) ? "Heatmap for sample pair: $sample_pair" : "Expression Heatmap for " . basename($sqlite_db);

    my $loader_name = "load_plots_$$";

    print $cgi->start_html(-title => $title,
                           -onLoad => $loader_name . "();",
                           -style => {'src' => "CSS/common.css"},
        );    


    print "<h1>$title</h1>\n";
    
    my $plot_loader = new CanvasXpress::PlotOnLoader($loader_name);

    my $heatmap_plot_obj = &write_heatmap($sqlite_db, \%params);

    if ($heatmap_plot_obj) {
        $plot_loader->add_plot($heatmap_plot_obj);
        
        print $plot_loader->write_plot_loader();
        
        print "<p style='clear:both;'>--</p>\n"; #just separate the annot divs that show up in the IGV go script.
        print &DEWebCommon::write_IGV_go_script($sqlite_db);
    }
    
    print $cgi->end_html();
    
    exit(0);
}

####
sub write_heatmap {
    my ($sqlite_db, $params_href) = @_;
    
    my $dbproc = connect_to_db($sqlite_db);
    

    ## none of this helps with performance...  :(
    #&AutoCommit($dbproc, 0);
    #&RunMod($dbproc, "PRAGMA synchronous=OFF");
    #&RunMod($dbproc, "pragma cache_size=4000000");
    #&RunMod($dbproc, "PRAGMA temp_store=MEMORY");
    #&RunMod($dbproc, "pragma journal_mode=memory");
    
    
    my $min_FC = $params_href->{min_FC} || 10;
    my $max_FDR = $params_href->{max_FDR} || 1e-10;

    my $scale_range = $params_href->{scale_range} || "min-max";

    


    if ($params_href->{all_features}) {
        $min_FC = undef;
        $max_FDR = undef;
    }
    
    my $min_any_feature_expr = $params_href->{min_any_feature_expr} || 0;
    my $min_sum_feature_expr = $params_href->{min_sum_feature_expr} || 0;
    

    my $sample_pair = $params_href->{sample_pair};
    my $restrict_to_sample_pair_flag = $params_href->{restrict_to_sample_pair_flag};


    my $max_genes_show = $params_href->{max_genes_show} || 100;    
    my $center_vals = $params_href->{center_vals} || "none";

    my $top_expressed_flag = $params_href->{top_expressed_flag} || 0;
    
    # print "<pre>\n" . Dumper($params_href) . "</pre>\n";
    
    my $html = "<div id='control_box'>\n";
    
    $html .= "<form id='control_form' method='get' action='HeatmapNav.cgi' >\n";
    $html .= "<ul>\n";
    $html .= "   <li>min_FC: <input type='text' name='min_FC' value=\"$min_FC\" />\n";
    $html .= "   <li>max_FDR: <input type='text' name='max_FDR' value=\"$max_FDR\" />\n";

    $html .= "   <li>min_any_expr_per_gene: <input type='text' name='min_any_feature_expr' value=\"$min_any_feature_expr\" />\n";
    $html .= "   <li>min_sum_feature_expr: <input type='text' name='min_sum_feature_expr' value=\"$min_sum_feature_expr\" />\n";
    
    $html .= "   <li>Heatmap scale range: <input type='text' name='scale_range' value=\'$scale_range\' size=8 />\n"; 
    
    $html .= "<li>Center expression values: ";

    my @center_val_opts = qw(average median none);
    foreach my $opt (@center_val_opts) {
        
        $html .= "<input type='radio' name=\'center_vals\' value=\'$opt\' ";
        if ($opt eq $center_vals) {
            $html .= " checked ";
        }
        $html .= " />$opt\n";
    }
    
    
    
    my $feature_type = $params_href->{feature_type} || "G"; # genes by default
    {
        # Feature type selection
    
        $html .= "<li>Feature type: ";
        $html .= "<input type='radio' name='feature_type' value='G' ";
        if ($feature_type eq 'G') {
            $html .= " checked ";
        }
        $html .= " />Genes\n";
        $html .= "<input type='radio' name='feature_type' value='T' ";
        if ($feature_type eq 'T') {
            $html .= " checked ";
        }
        $html .= " />Transcripts\n";
    }

    {

        $html .= "   <li><input type='checkbox' name='all_features' value='1' ";
        if ($params_href->{all_features}) {
            $html .= " checked ";
        }
        
        $html .= "/>All features (ignore min_FC, max_FDR)\n";
        
    }

    {
        $html .= "   <li><input type='checkbox' name='cluster_transcripts' value='1' ";
        if ($params_href->{cluster_transcripts}) {
            $html .= " checked ";
        }
        $html .= " />Cluster transcripts\n";
    }
    
    {
        $html .= "   <li><input type='checkbox' name='cluster_samples' value='1' ";
        if ($params_href->{cluster_samples}) {
            $html .= " checked ";
        }
        $html .= " />Cluster samples\n"; 
    }
    
    {
        $html .= "   <li><input type='checkbox' name='top_expressed_flag' value='1' ";
        if ($top_expressed_flag) {
            $html .= " checked ";
        }
        $html .= " />Restrict to top-most expressed in any given sample.\n";
    }
   
    if ($sample_pair) {
        $html .= "<input type='hidden' name='sample_pair' value=\'$sample_pair\' />\n";
        
        $html .= "<li><input type='checkbox' name='restrict_to_sample_pair_flag' value='1' ";
        if ($restrict_to_sample_pair_flag) {
            $html .= " checked ";
        }
        $html .= " />Restrict heatmap to sample pair\n";
                
    }
    
    $html .= "<li>Max genes to show: <input type='text' name='max_genes_show' value=$max_genes_show size=4 />\n";
       

    $html .= "   <li><input type='submit' />\n";
    $html .= "</ul>\n";

    $html .= "<input type='hidden' name='sqlite' value=\'$sqlite_db\' />\n";
    
    $html .= "</form>\n";
    $html .= "</div>\n";

    print $html;

    if ($center_vals && $center_vals eq 'none') {
        undef($center_vals);
    }    
    

    my %data = &DEWebCommon::get_expression_data($dbproc, { min_FC => $min_FC,
                                                            max_FDR => $max_FDR,
                                                            min_any_feature_expr => $min_any_feature_expr,
                                                            min_sum_feature_expr => $min_sum_feature_expr,
                                                            
                                                            sample_pair => $sample_pair, # just in case specified
                                                            restrict_to_sample_pair_flag => $restrict_to_sample_pair_flag,
                                                            
                                                            feature_type => $feature_type,
    
                                                            center_vals => $center_vals,
                                                            
                                                            max_select => 200000,
                                                 });
    
    #print Dumper(\%data);

    #######################
    # Downsample as needed
    #######################
    
        
    my @feature_names = keys %{$data{feature_to_info}};
    
    if (scalar @feature_names > $max_genes_show) {
        %data = &DEWebCommon::sample_from_data(\%data, $max_genes_show, $top_expressed_flag);
        
        print "<p> (Only $max_genes_show of " . scalar(@feature_names) . " randomly selected features are shown) </p>";
    }
    
    my $replicate_names_aref = $data{replicate_names};
    my $value_matrix_aref = $data{value_matrix};

    my $num_features = scalar(@$value_matrix_aref);
    
    print "<p>Found $num_features features.</p>\n";
    
    unless ($num_features) {
        # nothing to plot.
        return(undef);
    }
    
    
    my %heatmap_inputs = (

        sample_annotations => { samples => $data{sample_names} },
        
        replicate_names => $replicate_names_aref,
        value_matrix => $value_matrix_aref,
        cluster_features => $params_href->{cluster_transcripts},
        cluster_samples => $params_href->{cluster_samples},
        dendrogramSpace => 0.2,
        
                          
                          events => {

                            'click' => #"console.log('clicked' + o);\n"
                                 "var gene = o['y']['vars'][0];\n"
                                #. "console.log(gene);\n"
                                #. "alert('clicked ' + gene);\n"
                                . "IGV_go(gene);\n"
                                ,
                            #    
                            #'dblclick' => "console.log('dblclick' + Dumper.alert(o));\n",
                            #'dblclick' => "console.log(Dumper.write(cx.varDendrogram.nodes));\n",
#"Dumper.alert(o);\n",
                            'dblclick' => "var gene = o['y']['vars'][0];\n"
                                         . "launch_feature_report(gene);\n",

                            #
                            #    'mousemove' => #"DumperAlert(o);",
                                #"var gene = o['vars'][0];\n"
                            #    "console.log('1');\n"
                            #    . "alert('mouseover');\n",
                            
                          },
                          
        );
    

    # heatmap color scaling
    
    my ($min_scale, $max_scale) = split(/-/, $scale_range);
    if ($min_scale =~ /^\d+$/ && $max_scale =~ /^\d+$/) {
        $heatmap_inputs{setMinX} = $min_scale;
        $heatmap_inputs{setMaxX} = $max_scale;
    }
    
    my $heatmap_obj = new CanvasXpress::Heatmap("heatmap_$$");

    print $heatmap_obj->draw(%heatmap_inputs);

    
    return($heatmap_obj);
}

