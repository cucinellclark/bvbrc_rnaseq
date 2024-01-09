#
# The RNASeq Analysis application.
#

use strict;
use Carp;
use Data::Dumper;
use File::Temp;
use File::Slurp;
use File::Basename;
use IPC::Run 'run';
use File::Path 'rmtree';
use JSON;
use Bio::KBase::AppService::AppConfig;
use Bio::KBase::AppService::AppScript;
use Cwd;

our $global_ws;
our $global_token;

our $shock_cutoff = 10_000;

my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
# my $data_url = "http://www.alpha.patricbrc.org/api";

my $script = Bio::KBase::AppService::AppScript->new(\&process_rnaseq, \&preflight);
my $rc = $script->run(\@ARGV);
exit $rc;

# use JSON;
# my $temp_params = JSON::decode_json(`cat /home/fangfang/P3/dev_container/modules/app_service/test_data/rna.inp`);
# process_rnaseq('RNASeq', undef, undef, $temp_params);

# Return value: moved to global scope to work properly
my $run_ret_val = 0;

# flag for if localize_params is called, will remove the localized params after run finishes and before submitting to users workspace
my $called_localize_params = 0;

# flag for disabling uploading data to a user's workspace
my $disable_workspace_upload = 0;

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;
    my $mem_req = check_memory_requirements($app,$params);
    my $pf = {
	cpu => 8,
	memory => $mem_req,
	runtime => 0,
	storage => 0,
	is_control_task => 0,
    };
    return $pf;
}

#check memory requirements for preflight
#if above a certain threshold, set to 128GB
#otherwise, set memory requirements to 32GB
sub check_memory_requirements 
{
   my ($app,$params) = @_;
   my $mem_threshold = 50000000000; #50GB 
   my $total_mem = 0;
   my $ws = $app->workspace;
   my $comp_factor = 3; # compression factor
   #paired_end libs
   foreach my $item (@{$params->{paired_end_libs}}) {
      my $r1 = $ws->stat($item->{read1});
      my $r2 = $ws->stat($item->{read2});
      my $r1_size = $r1->size;
      my $r2_size = $r2->size;
      if ($ws->file_is_gzipped($item->{read1})) {
        $r1_size = $comp_factor*$r1_size;
      }
      if ($ws->file_is_gzipped($item->{read2})) {
        $r2_size = $comp_factor*$r2_size;
      }
      $total_mem = $total_mem + $r1_size + $r2_size;
   }
   #single_end libs
   foreach my $item (@{$params->{single_end_libs}}) { 
      my $r = $ws->stat($item->{read});
      my $r_size = $r->size;
      if ($ws->file_is_gzipped($item->{read})) {
        $r_size = $comp_factor*$r_size;
      }
      $total_mem = $total_mem + $r_size;
   }
   #bam libs
   foreach my $item (@{$params->{bam_libs}}) {
      my $b = $ws->stat($item->{bam});
      $total_mem = $total_mem + $b->size;
   }
   #check memory requirement and return 
   if ($total_mem >= $mem_threshold) {
      return "128GB";     
   } else {
      return "32GB";
   }
}

sub process_rnaseq {
    my ($app, $app_def, $raw_params, $params) = @_;
    print "Proc RNASeq ", Dumper($app_def, $raw_params, $params);
    my $time1 = `date`;

    my $parallel = $ENV{P3_ALLOCATED_CPU};

    #
    # Redirect tmp to large NFS if more than 4 input files.
    # (HACK)
    #
    my $file_count = count_params_files($params);
    print STDERR "File count: $file_count\n";
    my $bigtmp = "/vol/patric3/tmp";
    if ($file_count > 4 && -d $bigtmp)
    {
	print STDERR "Changing tmp from $ENV{TEMPDIR} to $bigtmp\n";
	$ENV{TEMPDIR} = $ENV{TMPDIR} = $bigtmp;
    }
    
    $global_token = $app->token();
    $global_ws = $app->workspace;
    
    my $output_folder = $app->result_folder();
    # my $output_base   = $params->{output_file};
    
    my $recipe = $params->{recipe};
    
    # my $tmpdir = File::Temp->newdir();
    my $tmpdir = File::Temp->newdir( CLEANUP => 1 );
    # my $tmpdir = File::Temp->newdir( CLEANUP => 0 );
    # my $tmpdir = "/tmp/RNApubref";
    # my $tmpdir = "/tmp/RNAuser";
    system("chmod", "755", "$tmpdir");
    print STDERR "$tmpdir\n";
    ###localize_params for regular script
    #localize_params_local for testing: will not download files
    $params = localize_params($tmpdir, $params);
    #$params = localize_params_local($tmpdir, $params);


    # 
    # If experimental_conditions was not specified, default to the empty list.
    #

    $params->{experimental_conditions} //= [];
    $params->{contrasts} //= [];
    
    my @outputs;
    my $prefix = $recipe;
    my $host = 0;
    if ($recipe eq 'cufflinks' || $recipe eq 'HTSeq-DESeq') {
        @outputs = run_bvbrc_rnaseq($params, $tmpdir, $host, $parallel);
    } elsif ($recipe eq 'Host') {
        $host = 1;
        @outputs = run_bvbrc_rnaseq($params, $tmpdir, $host, $parallel);
    } else {
        die "Unrecognized recipe: $recipe \n";
    }
    print STDERR '\@outputs = '. Dumper(\@outputs);

    # remove localized params (just the downloaded read files)
    if ($called_localize_params) {
        remove_localized_params($tmpdir, $params); 
    }

    if ($disable_workspace_upload) {
        die "disable_workspace_upload is true: terminating job before upload\n";
    }

    #
    # Create folders first.
    #
    for my $fent (grep { $_->[1] eq 'folder' } @outputs)
    {
	my $folder = $fent->[0];
	my $file = basename($folder);
	my $path = "$output_folder/$file";
	eval {
	    $app->workspace->create( { objects => [[$path, 'folder']] } );
	};
	if ($@)
	{
	    warn "error creating $path: $@";
	}
	else
	{
	    my $type ="txt";
	    if (opendir(my $dh, $folder))
	    {
		while (my $filename = readdir($dh))
		{
		    if ($filename =~ /\.json$/)
		    {
			my $ofile = "$folder/$filename";
			my $dest = "$path/$filename";
			print STDERR "Output folder = $folder\n";
			print STDERR "Saving $ofile => $dest ...\n";
			$app->workspace->save_file_to_file($ofile, {}, $dest, $type, 1,
							   (-s "$ofile" > $shock_cutoff ? 1 : 0), # use shock for larger files
							   $global_token);
		    }
		}
	    }
	    else
	    {
		warn "Cannot open output folder $folder: $!";
	    }
	}
    }
    for my $output (@outputs)
    {
	my($ofile, $type) = @$output;
	next if $type eq 'folder';
	
	if (! -f $ofile)
	{
	    warn "Output file '$ofile' of type '$type' does not exist\n";
	    next;
	}
	
	if ($type eq 'job_result')
	{
            my $filename = basename($ofile);
            print STDERR "Output folder = $output_folder\n";
            print STDERR "Saving $ofile => $output_folder/$filename ...\n";
            $app->workspace->save_file_to_file("$ofile", {},"$output_folder/$filename", $type, 1);
	}
    }
    my $time2 = `date`;
    my $outdir = "$tmpdir";
    
    if ($called_localize_params) {
        # Remove genome directory
        my $ref_id   = $params->{reference_genome_id};
        my $ref_dir  = "$tmpdir/$ref_id";
        if ( -d $ref_dir ) {
            rmtree([ "$ref_dir" ]);
        }

        # Remove TPMCalculator directory
        my $tpm_dir = "$tmpdir/TPMCalculator";
        if ( -d $tpm_dir ) {
            rmtree([ "$tpm_dir" ]);
        }
    }

    # save diffexp stuff
    #my $diffexp_name = "diff_exp";
    #my $diffexp_folder = "$output/.$diffexp_name";
    #my $diffexp_file = "$output/$diffexp_name";

    # copy over diffexp stuff
    #if ( -e $diffexp_folder ) {
    #    $app->workspace->save_file_to_file("$diffexp_folder", {},"$output_folder/.$diffexp_name", "job_result", 1);
    #}

    save_output_files($app,$outdir);
    write_output("Start: $time1"."End:   $time2", "$tmpdir/DONE");

    if (!$run_ret_val) {
	    die "Error running prok_tuxedo.py: saved any output to user's workspace\n";
    }
}

sub run_bvbrc_rnaseq {
    my ($params, $tmpdir, $host, $parallel) = @_;

    #$parallel //= 1;
    $parallel //= 8;
    
    my $cwd = getcwd();
    
    my $json = JSON::XS->new->pretty(1);

    #
    # Write job description.
    #


    my $jdesc = "$cwd/jobdesc.json";
    write_file($jdesc, $json->encode($params));
    
    my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
    my $dat = { data_api => "$data_api/genome_feature" };
    my $override = {
	cufflinks => { -p => $parallel},
	cuffdiff => {-p => $parallel},
	cuffmerge => {-p => $parallel},
	hisat2 => {-p => $parallel},
	bowtie2 => {-p => $parallel},
	stringtie => {-p => $parallel},
    htseq => {-p => $parallel}
    };
    #
    # no pretty, ensure it's on one line
    #i
    my $pstring = encode_json($override);
    my $sstring = encode_json($dat);
    
    my $outdir = "$tmpdir";
    
    my %exps     = params_to_exps($params);
    my $labels   = $params->{experimental_conditions};
    my $ref_id   = $params->{reference_genome_id} or die "Reference genome is required\n";
    my $output_name = $params->{output_file} or die "Output name is required\n";
    my $host_ftp = defined($params->{host_ftp}) ? $params->{host_ftp} : undef;
    my $dsuffix = "_diffexp";
    my $diffexp_name = ".$output_name$dsuffix";
    my $diffexp_folder = "$outdir/$diffexp_name";
    my $diffexp_file = "$outdir/$output_name$dsuffix";
    my $ref_dir  = prepare_ref_data_rocket($ref_id, $tmpdir, $host, $host_ftp);
    #my $unit_test = defined($params->{unit_test}) ? $params->{unit_test} : undef;
    
    print "Run rna_rocket ", Dumper(%exps, $labels, $tmpdir);
    
    # my $rocket = "/home/fangfang/programs/Prok-tuxedo/prok_tuxedo.py";
    my $rocket = "run_rnaseq";
    verify_cmd($rocket);
    
    my @cmd = ($rocket);
    #if ($host) {
    #    push @cmd, ("--index");
    #}
    push @cmd, ("-p", $pstring);
    push @cmd, ("-o", $outdir);
    push @cmd, ("-g", $ref_dir);
    push @cmd, ("-d", $diffexp_folder);
    push @cmd, ("--jfile", $jdesc);
    push @cmd, ("--sstring", $sstring);
    #if ($unit_test) {
    #    push @cmd, ("--unit_test",$params->{unit_test});
    #}
    
    #push @cmd, ("-L", join(",", map { s/^\W+//; s/\W+$//; s/\W+/_/g; $_ } @$labels)) if $labels && @$labels;
    #push @cmd, map { my @s = @$_; join(",", map { join("%", @$_) } @s) } @$exps;
    
    print STDERR "cmd = ", join(" ", @cmd) . "\n\n";
    
    #
    # Run directly with IPC::Run so that stdout/stderr can flow in realtime to the
    # output collection infrastructure.
    #
    #my $ok = run(\@cmd);
    $run_ret_val = run(\@cmd);
    if (!$run_ret_val)
    {
	    print "Error $? running @cmd\n";
    }
    
    #    my ($rc, $out, $err) = run_cmd(\@cmd);
    #    print STDERR "STDOUT:\n$out\n";
    #    print STDERR "STDERR:\n$err\n";
    
    run("echo $outdir && ls -ltr $outdir");
    
    #
    # Collect output and assign types.
    #
    my @outputs;

    # check for no_condition directory
    if (-d "$outdir/no_condition") {
        push @$labels, 'no_condition';
    }
    
    #
    # BAM/BAI/GTF files are in the replicate folders.
    # We flatten the file structure in replicate folders for the
    # files we are saving.
    #
    for my $cond (@$labels)
    {
	my @reps = map { basename($_) } glob("$outdir/$cond/*");
	
	for my $rep (@reps)
	{
	    my $path = "$outdir/$cond/$rep";
	    #
	    # Suffix/type list for output
	    #
	    my @types = (['.bam', 'bam'], ['.bai', 'bai'], ,['.gtf', 'gff'], ['.html', 'html'], ['_tracking', 'txt']);
	    for my $t (@types)
	    {
		my($suffix, $type) = @$t;
		for my $f (glob("$path/*$suffix"))
		{
		    my $base = basename($f);
			push(@outputs, [$f, $type]);
		}
	    }
	}
    }
    
    #
    # Remaining files are loaded as plain text.
    #
    for my $txt (glob("$outdir/*diff"))
    {
	push(@outputs, [$txt, 'txt']);
    }
    
    push @outputs, [ "$outdir/gene_exp.gmx", 'diffexp_input_data' ] if -s "$outdir/gene_exp.gmx";
    push @outputs, [ $diffexp_file, 'job_result' ] if -s $diffexp_file;
    push @outputs, [ $diffexp_folder, 'folder' ] if -e $diffexp_folder and -d $diffexp_folder;
    
    return @outputs;
}

sub run_rockhopper {
    my ($params, $tmpdir) = @_;
    
    my $exps     = params_to_exps($params);
    my $labels   = $params->{experimental_conditions};
    my $stranded = defined($params->{strand_specific}) && !$params->{strand_specific} ? 0 : 1;

    my $ref_id   = $params->{reference_genome_id};
    my $ref_dir  = prepare_ref_data($ref_id, $tmpdir) if $ref_id;
    
    print "Run rockhopper ", Dumper($exps, $labels, $tmpdir);
    
    # my $jar = "/home/fangfang/programs/Rockhopper.jar";
    my $jar = $ENV{KB_RUNTIME} . "/lib/Rockhopper.jar";
    -s $jar or die "Could not find Rockhopper: $jar\n";
    
    my $outdir = "$tmpdir/Rockhopper";
    
    my @cmd = (qw(java -Xmx1200m -cp), $jar, "Rockhopper");
    
    print STDERR '$exps = '. Dumper($exps);
    
    my @conditions = clean_labels($labels);
    
    push @cmd, qw(-SAM -TIME);
    push @cmd, qw(-s false) unless $stranded;
    push @cmd, ("-p", 1);
    push @cmd, ("-o", $outdir);
    push @cmd, ("-g", $ref_dir) if $ref_dir;
    push @cmd, ("-L", join(",", @conditions)) if $labels && @$labels;
    push @cmd, map { my @s = @$_; join(",", map { join("%", @$_) } @s) } @$exps;
    
    print STDERR "cmd = ", join(" ", @cmd) . "\n\n";
    
    my ($rc, $out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    
    run("echo $outdir && ls -ltr $outdir");
    
    my @outputs;
    if ($ref_id) {
        @outputs = merge_rockhoppper_results($outdir, $ref_id, $ref_dir);
        my $gmx = make_diff_exp_gene_matrix($outdir, $ref_id, \@conditions);
        push @outputs, [ $gmx, 'diffexp_input_data' ] if -s $gmx;
    } else {
        my @files = glob("$outdir/*.txt");
        @outputs = map { [ $_, 'txt' ] } @files;
    }
    
    return @outputs;
}

sub make_diff_exp_gene_matrix {
    my ($dir, $ref_id, $conditions) = @_;
    
    my $transcript = "$dir/$ref_id\_transcripts.txt";
    my $num = scalar@$conditions;
    return unless -s $transcript && $num > 1;
    
    my @genes;
    my %hash;
    my @comps;
    
    my @lines = `cat $transcript`;
    shift @lines;
    my $comps_built;
    for (@lines) {
        my @cols = split /\t/;
        my $gene = $cols[6]; next unless $gene =~ /\w/;
        my @exps = @cols[9..8+$num];
        # print join("\t", $gene, @exps) . "\n";
        push @genes, $gene;
        for (my $i = 0; $i < @exps; $i++) {
            for (my $j = $i+1; $j < @exps; $j++) {
                my $ratio = log_ratio($exps[$i], $exps[$j]);
                my $comp = comparison_name($conditions->[$i], $conditions->[$j]);
                $hash{$gene}->{$comp} = $ratio;
                push @comps, $comp unless $comps_built;
            }
        }
        $comps_built = 1;
    }
    
    my $outf = "$dir/$ref_id\_gene_exp.gmx";
    my @outlines;
    push @outlines, join("\t", 'Gene ID', @comps);
    for my $gene (@genes) {
        my $line = $gene;
        $line .= "\t".$hash{$gene}->{$_} for @comps;
        push @outlines, $line;
    }
    my $out = join("\n", @outlines)."\n";
    write_output($out, $outf);
    
    return $outf;
}

sub log_ratio {
    my ($exp1, $exp2) = @_;
    $exp1 = 0.01 if $exp1 < 0.01;
    $exp2 = 0.01 if $exp2 < 0.01;
    return sprintf("%.3f", log($exp2/$exp1) / log(2));
}

sub comparison_name {
    my ($cond1, $cond2) = @_;
    return join('|', $cond2, $cond1);
}

sub clean_labels {
    my ($labels) = @_;
    return undef unless $labels && @$labels;
    return map { s/^\W+//; s/\W+$//; s/\W+/_/g; $_ } @$labels;
}

sub merge_rockhoppper_results {
    my ($dir, $gid, $ref_dir_str) = @_;
    my @outputs;
    
    my @ref_dirs = split(/,/, $ref_dir_str);
    my @ctgs = map { s/.*\///; $_ } @ref_dirs;
    
    my %types = ( "transcripts.txt" => 'txt',
		 "operons.txt"     => 'txt' );
    
    for my $result (keys %types) {
        my $type = $types{$result};
        my $outf = join("_", "$dir/$gid", $result);
        my $out;
        for my $ctg (@ctgs) {
            my $f = join("_", "$dir/$ctg", $result);
            my @lines = `cat $f`;
            my $hdr = shift @lines;
            $out ||= join("\t", 'Contig', $hdr);
            $out  .= join('', map { join("\t", $ctg, $_ ) } grep { /\S/ } @lines);
        }
        write_output($out, $outf);
        push @outputs, [ $outf, $type ];
    }
    
    my @sams = glob("$dir/*.sam");
    for my $f (@sams) {
        my $sam = basename($f);
        my $bam = $sam;
        $bam =~ s/_R[12]\.sam$/.sam/;
        $bam =~ s/\.sam$/.bam/;
        $bam = "$dir/$bam";
        # my @cmd = ("samtools", "view", "-bS", $f, "-o", $bam);
        my @cmd = ("samtools", "sort", "-T", "$f.temp", "-O", "bam","-o", $bam, $f);
        run_cmd(\@cmd);
        push @outputs, [ $bam, 'bam' ];
        @cmd = ("samtools", "index", $bam);       
        run_cmd(\@cmd);
        push @outputs, [ "$bam.bai", 'bai' ];
    }
    push @outputs, ["$dir/summary.txt", 'txt'];
    
    return @outputs;
}

sub prepare_ref_data_rocket {
    my ($gid, $basedir, $host, $host_ftp) = @_;
    $gid or die "Missing reference genome id\n";
    
    my $dir = "$basedir/$gid";
    system("mkdir -p $dir");
    
    if ($host){
        if ($host_ftp){
            my $tar_url = "$host_ftp" . "_genomic.ht2.tar" ;
            my $out_file = basename $tar_url;
            my $out = curl_ftp($tar_url, "$dir/$out_file");
            my $fna_url = "$host_ftp" . "_genomic.fna" ;
            $out_file = basename $fna_url;
            $out = curl_ftp($fna_url, "$dir/$out_file");
            my $gff_url = "$host_ftp" . "_genomic.gff" ;
            $out_file = basename $gff_url;
            $out = curl_ftp($gff_url, "$dir/$out_file");
        }
        else{
            my $tar_url = "ftp://ftp.patricbrc.org/genomes/$gid/$gid.RefSeq.ht2.tar";
            my $out = curl_ftp($tar_url,"$dir/$gid.RefSeq.ht2.tar");
            my $fna_url = "ftp://ftp.patricbrc.org/genomes/$gid/$gid.RefSeq.fna";
            $out = curl_ftp($fna_url,"$dir/$gid.RefSeq.fna");
            my $gff_url = "ftp://ftp.patricbrc.org/genomes/$gid/$gid.RefSeq.gff";
            $out = curl_ftp($gff_url,"$dir/$gid.RefSeq.gff");
        }
    }
    
    else{
        my $api_url = "$data_url/genome_feature/?and(eq(genome_id,$gid),eq(annotation,PATRIC),eq(feature_type,CDS))&sort(+accession,+start,+end)&http_accept=application/cufflinks+gff&limit(25000)";
        #my $api_url = "$data_url/genome_feature/?and(eq(genome_id,$gid),eq(annotation,PATRIC),or(eq(feature_type,CDS),eq(feature_type,tRNA),eq(feature_type,rRNA)))&sort(+accession,+start,+end)&http_accept=application/cufflinks+gff&limit(25000)";
        my $ftp_url = "ftp://ftp.patricbrc.org/genomes/$gid/$gid.PATRIC.gff";
	
        my $url = $api_url;
        # my $url = $ftp_url;
        my $out = curl_text($url);
        write_output($out, "$dir/$gid.gff");

        # get list of valid accessions
        my @gff_lines = split(/\n/, $out);
        my %unique_accessions;
        foreach my $line (@gff_lines) {
            next if $line =~ /#/;
            my ($accession) = split(/\t/,$line);
            $unique_accessions{$accession} = 1;
        }
        my $accession_str = join(",", keys %unique_accessions);

        $api_url = "$data_url/genome_sequence/?eq(genome_id,$gid)&http_accept=application/sralign+dna+fasta&limit(25000)&in(accession,($accession_str))";
        $ftp_url = "ftp://ftp.patricbrc.org/genomes/$gid/$gid.fna";

        $url = $api_url;
        # $url = $ftp_url;
        my $out = curl_text($url);
        # $out = break_fasta_lines($out."\n");
        $out =~ s/\n+/\n/g;
        write_output($out, "$dir/$gid.fna");
    }
    
    return $dir;
}

sub prepare_ref_data {
    my ($gid, $basedir) = @_;
    $gid or die "Missing reference genome id\n";
    
    my $url = "$data_url/genome_sequence/?eq(genome_id,$gid)&select(accession,genome_name,description,length,sequence)&sort(+accession)&http_accept=application/json&limit(25000)";
    my $json = curl_json($url);
    # print STDERR '$json = '. Dumper($json);
    my @ctgs = map { $_->{accession} } @$json;
    my %hash = map { $_->{accession} => $_ } @$json;
    
    $url = "$data_url/genome_feature/?and(eq(genome_id,$gid),eq(annotation,PATRIC),eq(feature_type,CDS))&select(accession,start,end,strand,aa_length,patric_id,protein_id,gene,refseq_locus_tag,figfam_id,product)&sort(+accession,+start,+end)&limit(25000)&http_accept=application/json";
    $json = curl_json($url);
    
    for (@$json) {
        my $ctg = $_->{accession};
        push @{$hash{$ctg}->{cds}}, $_;
    }
    
    $url = "$data_url/genome_feature/?and(eq(genome_id,$gid),eq(annotation,PATRIC),or(eq(feature_type,tRNA),eq(feature_type,rRNA)))&select(accession,start,end,strand,na_length,patric_id,protein_id,gene,refseq_locus_tag,figfam_id,product)&sort(+accession,+start,+end)&limit(25000)&http_accept=application/json";
    $json = curl_json($url);
    
    for (@$json) {
        my $ctg = $_->{accession};
        push @{$hash{$ctg}->{rna}}, $_;
    }
    
    my @dirs;
    for my $ctg (@ctgs) {
        my $dir = "$basedir/$gid/$ctg";
        system("mkdir -p $dir");
        my $ent = $hash{$ctg};
        my $cds = $ent->{cds};
        my $rna = $ent->{rna};
        my $desc = $ent->{description} || join(" ", $ent->{genome_name}, $ent->{accession});
	
        # Rockhopper only parses FASTA header of the form: >xxx|xxx|xxx|xxx|ID|
        my $fna = join("\n", ">genome|$gid|accn|$ctg|   $desc   [$ent->{genome_name}]",
                       uc($ent->{sequence}) =~ m/.{1,60}/g)."\n";
	
        my $ptt = join("\n", "$desc - 1..$ent->{length}",
		       scalar@{$ent->{cds}}.' proteins',
		       join("\t", qw(Location Strand Length PID Gene Synonym Code FIGfam Product)),
		       map { join("\t", $_->{start}."..".$_->{end},
				  $_->{strand},
				  $_->{aa_length},
				  $_->{patric_id} || $_->{protein_id},
				  # $_->{refseq_locus_tag},
				  $_->{patric_id},
				  # $_->{gene},
				  join("/", $_->{refseq_locus_tag}, $_->{gene}),
				  '-',
				  $_->{figfam_id},
				  $_->{product})
			     } @$cds
                      )."\n" if $cds && @$cds;
	
        my $rnt = join("\n", "$desc - 1..$ent->{length}",
		       scalar@{$ent->{rna}}.' RNAs',
		       join("\t", qw(Location Strand Length PID Gene Synonym Code FIGfam Product)),
		       map { join("\t", $_->{start}."..".$_->{end},
				  $_->{strand},
				  $_->{na_length},
				  $_->{patric_id} || $_->{protein_id},
				  # $_->{refseq_locus_tag},
				  $_->{patric_id},
				  # $_->{gene},
				  join("/", $_->{refseq_locus_tag}, $_->{gene}),
				  '-',
				  $_->{figfam_id},
				  $_->{product})
			     } @$rna
                      )."\n" if $rna && @$rna;
	
        write_output($fna, "$dir/$ctg.fna");
        write_output($ptt, "$dir/$ctg.ptt") if $ptt;
        write_output($rnt, "$dir/$ctg.rnt") if $rnt;
	
        push(@dirs, $dir) if $ptt;
    }
    
    return join(",",@dirs);
}

sub curl_text {
    my ($url) = @_;
    my @cmd = ("curl", curl_options(), $url);
    print STDERR join(" ", @cmd)."\n";
    my ($out) = run_cmd(\@cmd);
    return $out;
}

sub curl_file {
    my ($url, $outfile) = @_;
    my @cmd = ("curl", curl_options(), "-o", $outfile, $url);
    print STDERR join(" ", @cmd)."\n";
    my ($out) = run_cmd(\@cmd);
    return $out;
}

sub curl_ftp {
    my ($url, $outfile) = @_;
    my @cmd = ("curl", "-o", $outfile, $url);
    print STDERR join(" ", @cmd)."\n";
    my ($out) = run_cmd(\@cmd);
    return $out;
}

sub curl_json {
    my ($url) = @_;
    my $out = curl_text($url);
    my $hash = JSON::decode_json($out);
    return $hash;
}

sub curl_options {
    my @opts;
    my $token = get_token()->token;
    push(@opts, "-H", "Authorization: $token");
    push(@opts, "-H", "Content-Type: multipart/form-data");
    return @opts;
}

sub run_cmd {
    my ($cmd) = @_;
    my ($out, $err);
    run($cmd, '>', \$out, '2>', \$err)
        or die "Error running cmd=@$cmd, stdout:\n$out\nstderr:\n$err\n";
    # print STDERR "STDOUT:\n$out\n";
    # print STDERR "STDERR:\n$err\n";
    return ($out, $err);
}

sub params_to_exps {
    my ($params) = @_;
    my %exps;
    for (@{$params->{paired_end_libs}}) {
        my $condition = $_->{condition};
        if (!exists($exps{$condition})) {
            $exps{$condition} = [];
        }
        push @{$exps{$condition}}, [ $_->{read1}, $_->{read2} ];
    }
    for (@{$params->{single_end_libs}}) {
        my $condition = $_->{condition};
        if (!exists($exps{$condition})) {
            $exps{$condition} = [];
        }
        push @{$exps{$condition}}, [ $_->{read} ];
    }
    return %exps;
}

# TODO: add removing srr libs
sub remove_localized_params {
    my ($tmpdir, $params) = @_;
    for (@{$params->{paired_end_libs}}) {
        unlink($_->{read1}) or warn "Can't delete $_->{read1}: $!\n";
        unlink($_->{read2}) or warn "Can't delete $_->{read2}: $!\n";
    }
    for (@{$params->{single_end_libs}}) {
        unlink($_->{read}) or warn "Can't delete $_->{read}: $!\n";
    }
}

sub localize_params {
    my ($tmpdir, $params) = @_;
    for (@{$params->{paired_end_libs}}) {
        $_->{read1} = get_ws_file($tmpdir, $_->{read1}) if $_->{read1};
        $_->{read2} = get_ws_file($tmpdir, $_->{read2}) if $_->{read2};
    }
    for (@{$params->{single_end_libs}}) {
        $_->{read} = get_ws_file($tmpdir, $_->{read}) if $_->{read};
    }
    $called_localize_params = 1;
    return $params;
}

sub localize_params_local {
    my ($tmpdir, $params) = @_;
    return $params;
}

sub count_params_files {
    my ($params) = @_;
    my $count = 0;
    if (ref($params->{paired_end_libs}))
    {
	$count += 2 * @{$params->{paired_end_libs}};
    }
    if (ref($params->{single_end_libs}))
    {
	$count += @{$params->{single_end_libs}};
    }
    return $count;
}

sub get_ws {
    return $global_ws;
}

sub get_token {
    return $global_token;
}

sub get_ws_file {
    my ($tmpdir, $id) = @_;
    # return $id; # DEBUG
    my $ws = get_ws();
    my $token = get_token();
    
    my $base = basename($id);
    my $file = "$tmpdir/$base";
    # return $file; # DEBUG
    
    my $fh;
    open($fh, ">", $file) or die "Cannot open $file for writing: $!";

    print STDERR "GET WS => $tmpdir $base $id\n";
    system("ls -la $tmpdir");

    eval {
	$ws->copy_files_to_handles(1, $token, [[$id, $fh]]);
    };
    if ($@)
    {
	die "ERROR getting file $id\n$@\n";
    }
    close($fh);
    print "$id $file:\n";
    system("ls -la $tmpdir");

    return $file;
}

sub write_output {
    my ($string, $ofile) = @_;
    open(F, ">$ofile") or die "Could not open $ofile";
    print F $string;
    close(F);
}

sub break_fasta_lines {
    my ($fasta) = @_;
    my @lines = split(/\n/, $fasta);
    my @fa;
    for (@lines) {
        if (/^>/) {
            push @fa, $_;
        } else {
            push @fa, /.{1,60}/g;
        }
    }
    return join("\n", @fa);
}

sub verify_cmd {
    my ($cmd) = @_;
    system("which $cmd >/dev/null") == 0 or die "Command not found: $cmd\n";
}

sub save_output_files
{
    my($app, $output) = @_;

    my $diffexp_name = "diff_exp";
    my $diffexp_folder = "$output/.$diffexp_name";
    my $diffexp_file = "$output/$diffexp_name";

    # copy over diffexp stuff
    #$app->workspace->save_file_to_file("$diffexp_folder", {},"$output/.$diffexp_name", "job_result", 1);
    
    my %suffix_map = (
              txt => 'txt',
              png => 'png',
              svg => 'svg',
              nwk => 'nwk',
              out => 'txt',
              err => 'txt',
              tsv => 'tsv',
              csv => 'csv',
              bam => 'bam',
              bai => 'bai',
              gtf => 'gff',
              gff => 'gff',
              _tracking => 'txt',
              gmx => 'diffexp_input_data',
              html => 'html');
           
    my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;

    if (opendir(my $dh, $output))
    {
    while (my $p = readdir($dh))
    {
        next if $p =~ /^\./;
        next if $p =~ /\.fna\z/;  
        next if $p eq $diffexp_folder;
        next if $p eq $diffexp_file;
        my @cmd = ("p3-cp", "-r", @suffix_map, "$output/$p", "ws:" . $app->result_folder);
        print "@cmd\n";
        my $ok = IPC::Run::run(\@cmd);
        if (!$ok)
        {
        warn "Error $? copying output with @cmd\n";
        }
    }
    closedir($dh);
    }
    else
    {
    warn "Output directory $output does not exist\n";
    }
}
