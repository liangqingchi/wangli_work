#!/usr/bin/perl

use File::Basename;

$num_args = $#ARGV + 1;

if ($num_args < 2) {
    print "\nUsage: \n  perl $0 <working_dir> <output_dir>\n  <working_dir>: store sample information, list and output.\n  <output_dir>: store the final 'clinical action' results in S3.\n\n";
    exit;
}
# find knowledge file which is in the same folder as script
my $dirname = dirname($0);

# find sample information file which is in working directory
my $workDir = $ARGV[0];
my $outDir = $ARGV[1];

use Data::UUID;

$ug = Data::UUID->new;
sub decrypt {
	my $encryptFile = shift @_;
	$uuid = $ug->create();
	$str = $ug->to_string($uuid);
	my $cmd = "echo ADMw3bmaster | gpg -o /var/tmp/$str.txt --passphrase-fd 0 $encryptFile";
	system($cmd);
	return "/var/tmp/$str.txt";
}


my $codeDescriptionPath = decrypt("$dirname/PGxOneV3_icd10_codes_desc.txt.gpg");
open (codeDescription, $codeDescriptionPath) || die "can't open ICD-9 Code Description file \n";


my %ICD_Code_Offical_Des;
my $line = <codeDescription>;

while ($line = <codeDescription>) {
	chomp($line);
	$line =~ s/\r//g;
	my @fields = split (/\t/, $line);
	$fields[0] =~ s/^\s//;
	$fields[0] =~ s/\s$//;
	$fields[1] =~ s/^\s//;
	$fields[1] =~ s/\s$//;
	$ICD_Code_Offical_Des{$fields[0]} = $fields[1];
}

close codeDescription;
system ("rm $codeDescriptionPath");

my $drugListPath = decrypt("$dirname/PGxOneV3_drug_list.txt.gpg");
open (drugList, $drugListPath) || die "can't open drug list file \n";

my %drug_weblink;

while ($line = <drugList>) {
	chomp($line);
	$line =~ s/\r//g;
	my @fields = split (/\t/, $line);
	$drug_weblink{lc $fields[0]} = $fields[1];
}
close drugList;
system ("rm $drugListPath");


my $drugActionPath = decrypt("$dirname/PGxOneV3_drug_action.txt.gpg");
open (drugAction, $drugActionPath) || die "can't open phenotype drug action file \n";

my $line = <drugAction>;
my %drug_phenotype;
my %drug_warning;
my %drug_recommendation;
my %drug_ICDCodes;
my %drug_Phenotype_description;
my %drug_link;
while ($line = <drugAction>) {
	chomp($line);
	$line =~ s/\r//g;
	my @fields = split (/\t/, $line);
	$fields[6] =~ s/^\s//;
	$fields[6] =~ s/\s$//;
	my $therapeutic_drug = $fields[0]."_".$fields[2]."_".$fields[3]."_".$fields[5];
	my $phenotype = $fields[5]."_".$fields[6];
	my $warning = $fields[1];
	my $recommendation = $fields[4];
	my $ICD_Codes = $fields[7];
	my $Phenotype_description = $fields[10];
	my $link = $fields[12];
	$drug_phenotype{$therapeutic_drug} = $phenotype;
	$drug_warning{$therapeutic_drug} =  $warning;
	$drug_recommendation{$therapeutic_drug} =  $recommendation;
	$drug_ICDCodes{$therapeutic_drug} = $ICD_Codes;
	$drug_Phenotype_description{$therapeutic_drug} = $Phenotype_description;
	$drug_link{$therapeutic_drug} = $link;
}
close drugAction;
system ("rm $drugActionPath");
open (sampleICDx, "$workDir/sample_codes_drugs.txt") || die "can't open sample ICD Codes and Rx drug list \n";
my %sample_ICDx;
my %sample_drugsTaken;

$line = <sampleICDx>;
while ($line = <sampleICDx>) {
	chomp($line);
	$line =~ s/\r//g;
	my @fields = split (/\t/, $line);
	my $sampleID = $fields[0];
	my $ICDx = $fields[1];
	$sample_ICDx{$sampleID} = $ICDx;
	$sample_drugsTaken{$sampleID} = $fields[2];
}
close sampleICDx; 	
open (sampleList, "$workDir/sample_variant_list.txt") || die "can't open sample list file \n";
$line = <sampleList>;
while ($line = <sampleList>) {
	chomp($line);
	$line =~ s/\r//g;
	my @fields = split (/\t/, $line);
	my $sampleID = $fields[1];
	my @temp1 = split (/_/,$sampleID);
	my $sampleID_singular = $temp1[0];
	if (($sampleID ne "PosCtrl-1") && ($sampleID ne "PosCtrl-2")) { 
		my $sample_phenotype_outputFile = $sampleID."_genotype_phenotype.txt";
		open (sample_phenotype, $sample_phenotype_outputFile) || die "can't open sample phenotype results file \n";
		
		my $output_file = ">$outDir/".$sampleID."_clinical_action.txt";
		open outputFile, $output_file;
		print outputFile "Therapeutic\tAction\tDrug Impacted\tClinical Interpretation\tGene\tGenotype\tPhenotype\tPhenotype_description\tlink\n";
		my %sample_gene_phenotype;
		my %sample_gene_phenotype_rs;  ###  update for seperate phenotype
		my %gene_phenotype;
		my %gene_phenotype_rs;  ###
		my %gene_genotype;
		my $CYP2C9_Code;
		my $genotype_phenotype_table = "";
		my $line1 = <sample_phenotype>;
		$genotype_phenotype_table = $genotype_phenotype_table.$line1;
		while ($line1 = <sample_phenotype>) {
			$genotype_phenotype_table = $genotype_phenotype_table.$line1;
			chomp($line1);
			$line1 =~ s/\r//g;
			my @fields1 = split (/\t/, $line1);
			my $gene = $fields1[0];
			my $genotype = $fields1[1];
			my $CYP2C9_note = $fields1[3];
			my $phenotype = $fields1[2];
			$phenotype =~ s/^\s//;
			$phenotype =~ s/\s$//;
			my @phenotype_set = split (/\//, $phenotype);
			my $count = scalar(@phenotype_set);
			for ($c=0;$c<$count;$c++) {  
				my $gene_phenotype =  $gene."_".$phenotype_set[$c];
				my  @phenotype_rs = split(/ /,$phenotype_set[$c]);                     ###
				my $gene_phenotype_rs =  $gene."_".$phenotype_rs[0];                   ###
				$sample_gene_phenotype{$gene_phenotype} = $genotype;
				$sample_gene_phenotype_rs{$gene_phenotype_rs} = $phenotype_set[$c];    ###
			
			}
			$gene_genotype{$gene} = $genotype;
			$gene_phenotype{$gene} = $phenotype;
			if ($gene eq "CYP2C9") {
				$CYP2C9_Code = $fields1[3];
			}
		}
		my $warfarin_rec, $warfarin_action;
		my $VOKRC1_genotype = $gene_genotype{"VKORC1"};
		if ($VOKRC1_genotype eq "WT/WT") {
			if (($CYP2C9_Code == 1) || ($CYP2C9_Code == 2)) {
				$warfarin_rec = "NORMAL DOSE Warfarin daily dose 5-7mg";
				$warfarin_action = "GO";
			}elsif (($CYP2C9_Code >= 3) || ($CYP2C9_Code <= 5)){
				$warfarin_rec = "DECREASE DOSE Warfarin daily dose 3-4mg";
				$warfarin_action = "DOWN";
			}else{
				$warfarin_rec = "DECREASE DOSE Warfarin daily dose 0.5-2mg";
				$warfarin_action = "DOWN";
			}
		}elsif ($VOKRC1_genotype eq "WT/-1639G>A") {
			if ($CYP2C9_Code == 1) {
				$warfarin_rec = "NORMAL DOSE Warfarin daily dose 5-7mg";
				$warfarin_action = "GO";
			}elsif (($CYP2C9_Code >= 2) || ($CYP2C9_Code <= 4)){
				$warfarin_rec = "DECREASE DOSE Warfarin daily dose 3-4mg";
				$warfarin_action = "DOWN";
			}else{
				$warfarin_rec = "DECREASE DOSE Warfarin daily dose 0.5-2mg";
				$warfarin_action = "DOWN";
			}
		}else {
			if (($CYP2C9_Code >= 1) || ($CYP2C9_Code <= 2)){
				$warfarin_rec = "DECREASE DOSE Warfarin daily dose 3-4mg";
				$warfarin_action = "DOWN";
			}else{
				$warfarin_rec = "DECREASE DOSE Warfarin daily dose 0.5-2mg";
				$warfarin_action = "DOWN";
			}
		}
				
		my %drug_action;
		my %current_drug_rec;
		my %current_drug_phe_des;
		my %current_drug_link;
		my %output_drug_phenotype;
		foreach my $drug (keys %drug_phenotype) {
			
			my @gene_phenotypes = split (/_/,$drug_phenotype{$drug});
			my $gene = $gene_phenotypes[0];
			my @phenotypes = split (/\//,$gene_phenotypes[1]);
			my @warnings = split (/\//,$drug_warning{$drug});
			my @drug_recs = split (/\//,$drug_recommendation{$drug});
			my @phe_dess = split (/\//,$drug_Phenotype_description{$drug});
			my $link = $drug_link{$drug};
			my $num_actionable_phenotypes = scalar(@phenotypes);
			my $action, $rec,$phe_des;
			my $flag = 0;
			my @drug_fields = split (/_/,$drug);
			my $drug_name = $drug_fields[2];
			if ($drug_name =~ m/^Warfarin/) {
				$action = $warfarin_action;
				$rec = $warfarin_rec;
				$flag =1;
			}else {
				for ($i = 0; $i<$num_actionable_phenotypes; $i++) {
					my $gene_phenotype = $gene."_".$phenotypes[$i];
					my @phenotype_rs = split(/ /,$phenotypes[$i]);    ###
					my $gene_phenotype_rs =  $gene."_".$phenotype_rs[0];   ###
					if (exists($sample_gene_phenotype{$gene_phenotype})) {
						$action = $warnings[$i];
						$output_drug_phenotype{$drug} = $phenotypes[$i];
						$phe_des = $phe_dess[$i];
						my @sub_recs = split (/;/,$drug_recs[$i]);
						$rec = $sub_recs[0];
						$flag =1;
						last;
					}
					if (not exists($sample_gene_phenotype{$gene_phenotype}) and exists($sample_gene_phenotype_rs{$gene_phenotype_rs})){ ###
						$output_drug_phenotype{$drug} = $sample_gene_phenotype_rs{$gene_phenotype_rs};                                  ###
					}
				}
			}
			if($flag == 0 ){
				$action = "GO";
				$rec = "NORMAL RESPONSE EXPECTED";
			}
			if ($action =~ m/^STOP/) {
				$action = "1_".$action;
			}
		
			if ($action =~ m/^DOWN/) {
				$action = "2_".$action;
			}
			if ($action =~ m/^UP/) {
				$action = "3_".$action;
			}
			if ($action =~ m/^SLOW/) {
				$action = "4_".$action;
			}
			if ($action =~ m/^GO/) {
				$action = "5_".$action;
			}
			my @drugs = split (/_/,$drug);
			my $therapeutic = $drugs[0];
			$action = $therapeutic."_".$action;
			
			$drug_action{$drug} = $action;
			$current_drug_rec{$drug} = $rec;
			$current_drug_phe_des{$drug} = $phe_des;
			$current_drug_link{$drug} = $link;
		}
		
		my %custom_drugs;
		my %sample_drugsTaken_rec;
		my $TCA_flag = 0;
		my %s_drugs;
		foreach my $drug (sort {$drug_action{$a} cmp $drug_action{$b} or lc $a cmp lc $b} keys %drug_action) {
			my @drugs = split (/_/,$drug);
			my $therapeutic = $drugs[0];
			my $drugImpacted = $drugs[1].":".$drugs[2];
			my $drug_name = $drugs[2];
			my @actions = split (/_/,$drug_action{$drug});
			my $action = $actions[2];
			my $rec = $current_drug_rec{$drug};
			my $phe_des = $current_drug_phe_des{$drug};
			my $link = $current_drug_link{$drug};
			my @gene_phenotype = split (/_/,$drug_phenotype{$drug});
			my $gene = $gene_phenotype[0];
			my $genotype = $gene_genotype{$gene};
			if(exists $output_drug_phenotype{$drug}){             ###
				$phenotype = $output_drug_phenotype{$drug};      ### 
			}                                                     ###
			else{
				@phenotype_wt = split(/\//,$gene_phenotype{$gene});
				$phenotype = $phenotype_wt[0];                   ###
			}                                                     ###
			if (($action eq "GO") &&($rec eq "NORMAL RESPONSE EXPECTED")){
				if (($phenotype =~ m/.*Metabolizer/) || ($phenotype =~ m/.*Expresser/)) {
					$phe_des = $phenotype
				}else{
					$phe_des = "Normal Response"
				}
			}
			if(($drug_name eq "Folic Acid") || ($drug_name eq "Pemetrexed (Alimta®)")){
				if($genotype !~m/.*C677T/) {
				$genotype = "WT/WT";
				$phenotype = "C677T Wild Type";
				}
			}
			print outputFile $therapeutic, "\t", $action, "\t", $drugImpacted, "\t", $rec, "\t", $gene, "\t", $genotype, "\t", $phenotype,"\t",$phe_des,"\t",$link,"\n";
			my @sample_ICD_Codes = split (/,/, (uc $sample_ICDx{$sampleID}));
			my $num_sample_codes = scalar(@sample_ICD_Codes);
			my %ICDx;

			for ($j=0;$j<$num_sample_codes; $j++){
					$sample_ICD_Codes[$j] =~ s/^\s//;
					$sample_ICD_Codes[$j] =~ s/\s$//;
				my @temp = split (/\./,$sample_ICD_Codes[$j]);
				my $digital_code = $temp[0];
				
				$ICDx{$digital_code} = 1;
				}
			
			my $ICD10 = $drug_ICDCodes{$drug};
			if ($ICD10 ne ""){
				my @codes = split (/,/,$ICD10);
				my $num_codes = scalar (@codes);
				my $code_flag = 0;
				for ($i=0;$i<$num_codes; $i++){
					$codes[$i] =~ s/^\s//;
					$codes[$i] =~ s/\s$//;
					if (exists $ICDx{$codes[$i]}){
						$code_flag = 1;
					}
				}
				if ($code_flag == 1) {
					$custom_drugs{$drug} = $drug_name;
				}
			}
			my @sample_drugs_taken = split (/,/, (lc $sample_drugsTaken{$sampleID}));
			my $num_sample_drugs = scalar(@sample_drugs_taken);
			for ($k=0;$k<$num_sample_drugs; $k++){
				$sample_drugs_taken[$k] =~ s/^\s//;
				$sample_drugs_taken[$k] =~ s/\s$//;
				$s_drugs{$sample_drugs_taken[$k]} = 1;
			}
			
			my @action_drugs = split (/,/, $drug_name);
			my $num_action_drugs = scalar (@action_drugs);
			my $drug_flag = 0;
			
			for ($p=0;$p<$num_action_drugs; $p++){
				my $length1 = index($action_drugs[$p], "(");
				my $length2 = index($action_drugs[$p], ")");
				my $generic_name = lc (substr $action_drugs[$p], 0, $length1-1);
				my $brand_name = lc (substr $action_drugs[$p], ($length1+1), ($length2-$length1-2)); 	
					$generic_name =~ s/^\s//;
					$generic_name =~ s/\s$//;
				$brand_name =~ s/^\s//;
				$brand_name =~ s/\s$//;
				my @generic_name_parts = split (/\//, $generic_name);
				my $num_generic_parts = scalar (@generic_name_parts);

					if (exists $s_drugs{$brand_name}){
						$drug_flag = 1;
					$sample_drugsTaken_rec{$action_drugs[$p]} =1;
					}
				for ($q=0;$q<$num_generic_parts; $q++){
					if (exists $s_drugs{$generic_name_parts[$q]}){
						$drug_flag = 1;
						$sample_drugsTaken_rec{$action_drugs[$p]} = 1;
					}
				}
			}
			if ($drug_flag == 1) {	
				if(!exists $custom_drugs{$drug}){
					$custom_drugs{$drug} = $drug_name;
				}
				}
		}
		print outputFile "\n";
		print outputFile "ICD-10:";
		my @sample_ICD_Codes = split (/,/,(uc $sample_ICDx{$sampleID}));
		my $num_sample_codes = scalar(@sample_ICD_Codes);
			for ($j=0;$j<$num_sample_codes; $j++){
				$sample_ICD_Codes[$j] =~ s/^\s//;
				$sample_ICD_Codes[$j] =~ s/\s$//;
			my @temp = split(/\./, $sample_ICD_Codes[$j]);
				my $full_code = $temp[0].$temp[1];
			if ($full_code ne "NA") {
				print outputFile " ",$sample_ICD_Codes[$j], " \"",$ICD_Code_Offical_Des{$full_code}, "\"";
			}
		}

		print outputFile "\n\n";
		print outputFile $genotype_phenotype_table;
		print outputFile "\n";

		my %stop_drugs;
		my $stop_drugs_rec = "";
		my %change_drugs;
		my %caution_drugs;
		my $caution_drugs_rec = "";
		my %go_drugs;
		my $go_drugs_rec = "";
		
		my @sample_drug_list = split(/,/, $sample_drugsTaken{$sampleID_singular});
		my %sample_drug_list_hash = map {$_ => 1} @sample_drug_list;
		
		foreach my $drug (sort {$custom_drugs{$a} cmp $custom_drugs{$b}} keys %custom_drugs) {
			my @drugs = split (/_/,$drug);
			my $drug_name = $drugs[2];			
			my @actions = split (/_/,$drug_action{$drug});
			my $action = $actions[2];
			my $rec = $current_drug_rec{$drug};

			if ($action =~ m/.*STOP.*/) {
				$stop_drugs{$drug_name} = 1;
				$stop_drugs_rec = "CONSIDER ALTERNATIVES";
			}
			if ($action =~ m/.*UP.*/) {
				my $used_rec;
				my $position = index($rec, "INCREASE DOSE");
				if ($position != -1) {
					my $new_rec = substr $rec, $position;
					
						if ((index($new_rec, "INCREASE DOSE by")) >=0) {
						my $end_pos = index($new_rec, "%") + 1;
							$used_rec = substr $new_rec, 0, $end_pos;
						}else{
							$used_rec = substr $new_rec, 0, 13;
					}
				}
				$change_drugs{$drug_name} = $used_rec;
			}

			if ($action =~ m/.*DOWN.*/) {
				my $used_rec;
				my $end_pos;
				my $position = index($rec, "DECREASE DOSE");
				if ($position != -1) {
					my $new_rec = substr $rec, $position;

					if ((index($new_rec, "DECREASE DOSE by")) >=0) {
						$end_pos = index($new_rec, "%") + 1;
					}elsif ((index($new_rec, "DECREASE DOSE to")) >=0){
						$end_pos = index($new_rec, "daily") + 5;
					}elsif ((index($new_rec, "DECREASE DOSE Warfarin")) >=0){
						$end_pos = index($new_rec, "mg") + 2;
					}else{
						$end_pos =13;
					}
						$used_rec = substr $new_rec, 0, $end_pos;
				}
					$change_drugs{$drug_name} = $used_rec;

			}

			if ($action =~ m/.*SLOW.*/) {
					$caution_drugs{$drug_name} = 1;
				$caution_drugs_rec = "USE CAUTION";
			}
				
			if ($action =~ m/.*GO.*/) {
				$go_drugs{$drug_name} = 1;
				$go_drugs_rec = "NORMAL RESPONSE EXPECTED";
			}
		}
		my $stop_drugs = "";
		if (!keys %stop_drugs){
			$stop_drugs = " ";
			$stop_drugs_rec = " ";
		}else{
			foreach my $ind_drug (sort keys %stop_drugs){
				$stop_drugs = $stop_drugs.$ind_drug.",";
			}
		}

		my $change_drugs = "";
		my $change_drugs_rec = "";
		if (!keys %change_drugs){
			$change_drugs = " ";
			$change_drugs_rec = " ";
		}else{
			foreach my $ind_drug (sort keys %change_drugs){
				$change_drugs = $change_drugs.$ind_drug.",";
				$change_drugs_rec = $change_drugs_rec.$change_drugs{$ind_drug}.",";
			}
		}

		my $go_drugs = "";
		if (!keys %go_drugs){
			$go_drugs = " ";
			$go_drugs_rec = " ";
		}else{
			foreach my $ind_drug (sort keys %go_drugs){
				$go_drugs = $go_drugs.$ind_drug.",";
			}
		}

		my $caution_drugs = "";
		if(!keys %caution_drugs){
			$caution_drugs  = " ";
			$caution_drugs_rec = " ";
		}else{
			foreach my $ind_drug (sort keys %caution_drugs){
				$caution_drugs = $caution_drugs.$ind_drug.",";
			}
		}

		print outputFile "STOP", "\t", $stop_drugs, "\t", $stop_drugs_rec, "\n";
		print outputFile "UP_DOWN", "\t", $change_drugs, "\t", $change_drugs_rec, "\n";
		print outputFile "CAUTION", "\t", $caution_drugs, "\t", $caution_drugs_rec,  "\n";
		print outputFile "GO", "\t", $go_drugs, "\t", $go_drugs_rec, "\n\n";

		print outputFile "DrugTaken:";
		foreach my $drug_taken (keys %sample_drugsTaken_rec){
			print outputFile "\t", $drug_taken;
		}
		print outputFile "\n";

		print outputFile "ZNA:";
		foreach my $temp_drug (keys %s_drugs) {
			print outputFile "\t", $temp_drug;
		}
		print outputFile "\n";
		my %drug_custom_sort_action;
		foreach my $drug (keys %drug_action) {
			if (exists $custom_drugs{$drug}){
				my @actions = split (/_/,$drug_action{$drug});
				my $sort_action;
				if ($actions[2] eq "STOP") {
					$sort_action = "01_STOP";
				}
				if ($actions[2] eq "STOP UP") {
					$sort_action = "02_STOP UP";
				}
				if ($actions[2] eq "STOP DOWN") {
					$sort_action = "03_STOP DOWN";
				}
				if ($actions[2] eq "STOP SLOW") {
					$sort_action = "04_STOP SLOW";
				}
				if ($actions[2] eq "UP") {
					$sort_action = "05_UP";
				}
				 if ($actions[2] eq "UP SLOW") {
					$sort_action = "06_UP SLOW";
				}
				if ($actions[2] eq "DOWN") {
					$sort_action = "07_DOWN";
				}
				if ($actions[2] eq "DOWN SLOW") {
					$sort_action = "08_DOWN SLOW";
				}
				if ($actions[2] eq "SLOW") {
					$sort_action = "09_SLOW";
				}
				if ($actions[2] eq "GO") {
					$sort_action = "10_GO";
				}
				my @drugs = split (/_/,$drug);
				my $sort_drug = $drugs[1]."_".$drugs[2]."_".$drugs[3]."_".$drugs[0]; 
				$drug_custom_sort_action{$sort_drug} = $sort_action;
			}
		}
		my $Atorvastatin_flag = 0;
		my $Lovastatin_flag = 0;
		my $Simvastatin_flag = 0;
		foreach my $sort_drug (sort {$drug_custom_sort_action{$a} cmp $drug_custom_sort_action{$b} or lc $a cmp lc $b} keys %drug_custom_sort_action) {
			my @drugs = split (/_/,$sort_drug);
			my $therapeutic = $drugs[3];
			my $drugImpacted = $drugs[0].":".$drugs[1];
			my $drug_name = $drugs[1];
			my @drug_set = split (/, /,$drug_name);
			my $ori_drug = $drugs[3]."_".$drugs[0]."_".$drugs[1]."_".$drugs[2];
			my @actions = split (/_/,$drug_custom_sort_action{$sort_drug});
			my $action = $actions[1];
			my $rec = $current_drug_rec{$ori_drug};
			my $phe_des = $current_drug_phe_des{$ori_drug};
			my $link = $current_drug_link{$ori_drug};
			my @gene_phenotype = split (/_/,$drug_phenotype{$ori_drug});
			my $gene = $gene_phenotype[0];
			my $genotype = $gene_genotype{$gene};
			my $phenotype = $gene_phenotype{$gene};
			my $output_flag = 1;
			if (($drug_name=~m/Simvastatin/) && ($Simvastatin_flag == 1)){
				pop (@drug_set);
			}
			if (($drug_name=~m/Atorvastatin.*/) && ($Atorvastatin_flag == 1)) {
				if(($drug_name=~m/Lovastatin/) && ($Lovastatin_flag == 1)){
					$output_flag = 0;
				}else{
					shift (@drug_set);
				}
			}
			my $drug_number = scalar(@drug_set);
			if (($output_flag==1) && ($drug_number!=0)){
				my $drugImpacted_original = $drugImpacted;
				my $new_drug_name = join(', ', @drug_set);
				$new_drug_name =~ s/^\s//;
				$new_drug_name =~ s/\s$//;
				$drugImpacted = $drugs[0].":".$new_drug_name;
				if (($action eq "GO") &&($rec eq "NORMAL RESPONSE EXPECTED")){
					if (($phenotype =~ m/.*Metabolizer/) || ($phenotype =~ m/.*Expresser/)){
						$phe_des = $phenotype
					}else{
						$phe_des = "Normal Response"
					}
				}
				print outputFile $therapeutic, "\t", $action, "\t", $drugImpacted, "\t", $rec, "\t", $gene, "\t", $genotype, "\t", $phenotype,"\t",$phe_des,"\t",$link,"\n";
				if ($action eq "GO") {
					my @all_drugs = split (/,/,$drug_name);
					my $drug_count = scalar(@all_drugs);
					my @web_links;
					for ($n=0;$n<$drug_count;$n++){
						my @temp = split (/\(/,(lc $all_drugs[$n]));
							$temp[0] =~ s/^\s//;
							$temp[0] =~ s/\s$//;
						my $generic_drug = $temp[0];
						if ((exists $drug_weblink{$generic_drug})&&($drug_weblink{$generic_drug} ne "")) {
							my $drug_website = $all_drugs[$n].": ".$drug_weblink{$generic_drug};
							push (@web_links, $drug_website);
						}
					}
					my $weblink_count = scalar(@web_links);
					for ($p=0;$p<$weblink_count-1;$p++){
#						print outputFile $web_links[$p], ",";
					}
#					print outputFile $web_links[$weblink_count-1],"\t";
				}else{
#					print outputFile "\t";
				}
#				print outputFile $drugImpacted_original, "\n";
			}
			if ($drug_name=~m/Atorvastatin/){
				$Atorvastatin_flag = 1;
			}

			if ($drug_name=~m/Lovastatin/){
				$Lovastatin_flag = 1;
			}
			if ($drug_name=~m/Simvastatin/){
				$Simvastatin_flag = 1;
			}
		}
	}
	close outputFile;
	close sample_phenotype;
	
}
close sampleList;

exit(0);
