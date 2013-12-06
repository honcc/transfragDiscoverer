#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Math::Random;
use Math::Round;
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a perl script to discover transcribed fragments by using optimized combination of gap size, split window size and split window fold.
#
#	Input
#
#		--queryCovPlsIndexPath=				path[compulsory]; path of the query index.hsh.pls that specify the perl pileups of read coverage; should be generated from pileupPerlStorableAntisenseArtefactCorrector;
#		--gffPath=							path[compulsory]; path of the reference GFF for gene annotation;
#		--outDir=							output directory
#
#	v0.1-0.5
#	[OBSOLETE]
#
#	v0.6
#	[Thu  5 Sep 2013 21:05:40 CEST] debut, based on cisNATFinder_v0.1;
#
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-09-19 12:30]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/transfragDiscoverer/v0.6/transfragDiscoverer_v0.6.pl --queryCovPlsIndexPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.95/correctedCovPls.upr/index.hsh.pls --gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff --outDir=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.95/transfragDiscoverer/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/transfragDiscoverer/v0.6/transfragDiscoverer_v0.6.pl
#	--queryCovPlsIndexPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.999/correctedCovPls.upr/index.hsh.pls
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff
#	--outDir=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.999/transfragDiscoverer/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $globalScriptDirPath = dirname(rel2abs($0));
open DEBUGLOG, ">", "$globalScriptDirPath/debug.log.txt";
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|1518, readParameters|1719
#	secondaryDependOnSub: currentTime|818
#
#<section ID="startingTasks" num="0">
#----------Read parameters ----------#
&printCMDLogOrFinishMessage("CMDLog");#->1518
my ($queryCovPlsIndexPath, $gffPath, $outDir) = &readParameters();#->1719
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedPath" num="1">
my $basicIGVTrackHsh_ref = {
	'GCPath' => "/Volumes/A_MPro2TB/NGS/IGVInfo/infoTracks/EHI/GC.Repetitiveness/win.100_step.10.GC.tdf", 
	'reptvPath' => "/Volumes/A_MPro2TB/NGS/IGVInfo/infoTracks/EHI/GC.Repetitiveness/win.100_step.10.Reptv.tdf", 
	'plusCovTDFPath' => "/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.95/correctedWig.upr/corrected.plus.wig.tdf",
	'minusCovTDFPath' => "/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.95/correctedWig.upr/corrected.minus.wig.tdf",
	'gffPath' => $gffPath,
};
my $IGVGenomeID = "EHI_v13";
my $IGVGffTrackHsh_ref = {};
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="2">
my $maxThread = 12;

#my $maxTransfragGapTorleranceAry_ref = [10, 15, 20, 25, 30, 35, 40, 45, 50];
#my $splitWinSizeAry_ref = [150, 125, 100, 75, 50];
#my $splitFoldAry_ref = [50, 40, 30, 20, 10];
#my $minTransfragLenAry_ref = [200, 150, 100];
#my $minTransfragPosCovAry_ref = [3, 2, 1];
#my $minTransfragCovPerNtAry_ref = [3, 2, 1];

my $maxTransfragGapTorleranceAry_ref = [20, 30, 40];
my $splitWinSizeAry_ref = [150, 125, 100, 75, 50];
my $splitFoldAry_ref = [50, 20];
my $minTransfragLenAry_ref = [150];
my $minTransfragPosCovAry_ref = [1];
my $minTransfragCovPerNtAry_ref = [2];
my $weightFactorHsh_ref = {

	#---negative factors, the less the better
	'NATTrnsfrgNum_AS'	=> -1, #---num of transfrg per mRNA with antisense transfrag
	'mRNATrnsfrgNum_SS'	=> -3, #---num of transfrg per mRNA with sense transfrag
	'neighborOvrlp_SS'	=> -3, #---percentage of mRNA with sense transfrag overlapping with it's neighbor
	'neighborOvrlp_AS'	=> -1, #---percentage of mRNA with antisense transfrag overlapping with it's neighbor

	#---postive factors, the more the better
	'pctWithTrnsfrg_AS'	=> 4, #---percentage of mRNA that were found to have antisense transfrags
	'pctWithTrnsfrg_SS'	=> 0, #---[optional], percentage of mRNA that were found to have sense transfrags
	'avgTrnfrgLen_AS'	=> 4, #---average length of sense transfrag
	'avgTrnfrgLen_SS'	=> 0, #---average length of sense transfrag
	'mRNACoverPct_SS'	=> 0, #---[optional], percentage of mRNA covered by sense tranfrag
};
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="3">
my @mkDirAry;
my $resultDir = $outDir;
my $resultGffDir = "$resultDir/gff/"; push @mkDirAry, $resultGffDir;
my $resultLogDir = "$resultDir/log/"; push @mkDirAry, $resultLogDir;
my $resultWigDir = "$resultDir/wig/"; push @mkDirAry, $resultWigDir;
my $resultXMLDir = "$resultDir/xml/"; push @mkDirAry, $resultXMLDir;
my $resultStorableDir = "$resultDir/storable/"; push @mkDirAry, $resultStorableDir;
my $ggplotDirHsh_ref = {};
foreach (qw /dat pdf R log/) {$ggplotDirHsh_ref->{$_} = "$resultDir/ggplot/$_"; push @mkDirAry, $ggplotDirHsh_ref->{$_};}
system ("mkdir -pm 777 $_") foreach @mkDirAry;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_defineOutFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutFilePath" num="4">
my $XMLPath = "$resultXMLDir/transfragDiscoverer.igv.xml";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_processInputData
#	primaryDependOnSub: checkGeneInfo|372, checkmRNAProximity|633, classifymRNABasedOnProximity|770, getCtgryGeneInfo|915, getIndivCntgCovPlsPath|949, getmRNASubset|1008, readGFF_oneRNAPerGene|1620
#	secondaryDependOnSub: checkOverlapAndProximity|398, currentTime|818, reportStatus|1752
#
#<section ID="processInputData" num="5">
#----------Read individual cntg pls
my $queryCovPlsPathHsh_ref = &getIndivCntgCovPlsPath($queryCovPlsIndexPath);#->949

#----------Read Gff
my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);#->1620
&checkGeneInfo($geneInfoHsh_ref);#->372

#----------Get bfmRNA and mRNA ranges
my @mRNAAry = qw/bfmRNA/;
my $mRNAInfoHsh_ref = &getCtgryGeneInfo($geneInfoHsh_ref, \@mRNAAry);#->915
my ($hitAndPrxmtyBymRNAHsh_ref) = &checkmRNAProximity($resultStorableDir, $geneInfoHsh_ref, $mRNAInfoHsh_ref, $maxThread);#->633
my ($mRNAPrxmtyClassListHsh_ref) = &classifymRNABasedOnProximity($hitAndPrxmtyBymRNAHsh_ref, $mRNAInfoHsh_ref, undef, undef);#->770
my ($mRNAInfoSubSetHsh_ref, $queryCovSubsetPlsPathHsh_ref) = &getmRNASubset($mRNAInfoHsh_ref, $mRNAPrxmtyClassListHsh_ref, $queryCovPlsPathHsh_ref);#->1008
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_optimizeParameters
#	primaryDependOnSub: calculateOptimalTransfragScore|276, optimizeTransfrgIdentificationParameters|1254, plotTransfrgIdentificationParameters|1424
#	secondaryDependOnSub: checkmRNATrnsfrgOverlap|664, evaluatemRNATrnsfrgOvrlp|837, ggplotXYLineThreeFactors|1039, identifyTransfrag|1078, printBothFHAndStdout|1495, reportStatus|1752
#
#<section ID="optimizeParameters" num="6">
my ($optTransfrgParamDataHsh_ref) = &optimizeTransfrgIdentificationParameters($queryCovSubsetPlsPathHsh_ref, $mRNAInfoSubSetHsh_ref, $maxTransfragGapTorleranceAry_ref, $splitWinSizeAry_ref, $splitFoldAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $minTransfragCovPerNtAry_ref, $maxThread, $resultStorableDir);#->1254
my ($optDataForPrintHsh_ref, $optDataForPlotHsh_ref, $optimalParamAry_ref, $dataTypeStatHsh_ref) = &calculateOptimalTransfragScore($optTransfrgParamDataHsh_ref, $weightFactorHsh_ref, $maxTransfragGapTorleranceAry_ref, $splitWinSizeAry_ref, $splitFoldAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $minTransfragCovPerNtAry_ref, $resultLogDir);#->276
&plotTransfrgIdentificationParameters($optDataForPlotHsh_ref, $ggplotDirHsh_ref, $optimalParamAry_ref, $maxTransfragGapTorleranceAry_ref, $splitWinSizeAry_ref, $splitFoldAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $minTransfragCovPerNtAry_ref, $dataTypeStatHsh_ref);#->1424
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_identifyTransfragsUsingOptParam
#	primaryDependOnSub: getTransfragUsingOptimalParameters|983
#	secondaryDependOnSub: identifyTransfrag|1078
#
#<section ID="identifyTransfragsUsingOptParam" num="7">
my ($trnsfrgInfoHsh_ref) = &getTransfragUsingOptimalParameters($optimalParamAry_ref, $queryCovPlsPathHsh_ref, $maxThread, $resultStorableDir);#->983
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 8_outputFiles
#	primaryDependOnSub: outputFinalTransfragGff|1311, outputIGVXML|1343
#	secondaryDependOnSub: printGFF_oneRNAPerGene_chooseStrnd_filter|1551
#
#<section ID="outputFiles" num="8">
&outputFinalTransfragGff($resultGffDir, $trnsfrgInfoHsh_ref, $IGVGffTrackHsh_ref);#->1311
&outputIGVXML($XMLPath, $basicIGVTrackHsh_ref, $IGVGenomeID, $IGVGffTrackHsh_ref);#->1343
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 9_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|1518
#	secondaryDependOnSub: currentTime|818
#
#<section ID="finishingTasks" num="9">
close DEBUGLOG;
&printCMDLogOrFinishMessage("finishMessage");#->1518
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	XML [n=1]:
#		outputIGVXML
#
#	general [n=8]:
#		checkGeneInfo, currentTime, getCtgryGeneInfo
#		printBothFHAndStdout, printCMDLogOrFinishMessage, readGFF_oneRNAPerGene
#		readParameters, reportStatus
#
#	gff [n=3]:
#		checkGeneInfo, printGFF_oneRNAPerGene_chooseStrnd_filter, readGFF_oneRNAPerGene
#
#	ggplot [n=1]:
#		ggplotXYLineThreeFactors
#
#	multithread [n=1]:
#		generateThreadHshWithRandomCntg
#
#	plotInR [n=1]:
#		ggplotXYLineThreeFactors
#
#	range [n=4]:
#		checkOverlapAndProximity, checkmRNAProximity, checkmRNATrnsfrgOverlap
#		classifymRNABasedOnProximity
#
#	reporting [n=2]:
#		currentTime, printBothFHAndStdout
#
#	specific [n=1]:
#		classifymRNABasedOnProximity
#
#	storable [n=1]:
#		getIndivCntgCovPlsPath
#
#	unassigned [n=9]:
#		calculateOptimalTransfragScore, evaluatemRNATrnsfrgOvrlp, getTransfragUsingOptimalParameters
#		getmRNASubset, identifyTransfrag, optimizeTransfrgIdentificationParameters
#		outputFinalTransfragGff, plotTransfrgIdentificationParameters, scanDramaticCoverageChange
#
#====================================================================================================================================================#

sub calculateOptimalTransfragScore {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: printBothFHAndStdout|1495
#	appearInSub: >none
#	primaryAppearInSection: 6_optimizeParameters|185
#	secondaryAppearInSection: >none
#	input: $maxTransfragGapTorleranceAry_ref, $minTransfragCovPerNtAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $optTransfrgParamDataHsh_ref, $resultLogDir, $splitFoldAry_ref, $splitWinSizeAry_ref, $weightFactorHsh_ref
#	output: $dataTypeStatHsh_ref, $optDataForPlotHsh_ref, $optDataForPrintHsh_ref, $optimalParamAry_ref
#	toCall: my ($optDataForPrintHsh_ref, $optDataForPlotHsh_ref, $optimalParamAry_ref, $dataTypeStatHsh_ref) = &calculateOptimalTransfragScore($optTransfrgParamDataHsh_ref, $weightFactorHsh_ref, $maxTransfragGapTorleranceAry_ref, $splitWinSizeAry_ref, $splitFoldAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $minTransfragCovPerNtAry_ref, $resultLogDir);
#	calledInLine: 191
#....................................................................................................................................................#
	my ($optTransfrgParamDataHsh_ref, $weightFactorHsh_ref, $maxTransfragGapTorleranceAry_ref, $splitWinSizeAry_ref, $splitFoldAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $minTransfragCovPerNtAry_ref, $resultLogDir) = @_;

	#---transform the optTransfrgParamDataHsh
	my $optDataForPrintHsh_ref = {};
	my $dataForScalingHsh_ref = {};
	my $dataTypeStatHsh_ref = {};
	my $optDataForPlotHsh_ref = {};

	foreach my $dataType (keys %{$optTransfrgParamDataHsh_ref}) {
		foreach my $minTransfragLen (sort {$a <=> $b} @{$minTransfragLenAry_ref}) {
			foreach my $minTransfragPosCov (sort {$a <=> $b} @{$minTransfragPosCovAry_ref}) {
				foreach my $minTransfragCovPerNt (sort {$a <=> $b} @{$minTransfragCovPerNtAry_ref}) {
					foreach my $maxTransfragGapTorlerance (sort {$a <=> $b} @{$maxTransfragGapTorleranceAry_ref}) {
						foreach my $splitWinSize (sort {$a <=> $b} @{$splitWinSizeAry_ref}) {
							foreach my $splitFold (sort {$a <=> $b} @{$splitFoldAry_ref}) {
								my $value = $optTransfrgParamDataHsh_ref->{$dataType}{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold};
								$optDataForPrintHsh_ref->{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold}{$dataType} = $value;
								push @{$dataForScalingHsh_ref->{$dataType}}, $value;
								$optDataForPlotHsh_ref->{$dataType}{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold} = $value;
							}
						}
					}
				}
			}
		}
	}
	
	#---calculate the mean for scaling
	foreach my $dataType (keys %{$dataForScalingHsh_ref}) {
		my $valueStatObj = Statistics::Descriptive::Full->new();
		$valueStatObj->add_data(@{$dataForScalingHsh_ref->{$dataType}});
		$dataTypeStatHsh_ref->{$dataType}{'mean'} = $valueStatObj->mean();
	}

	#---calculate the score
	#---the score = sum of the weighted and scaled deviation from mean
	my $optimalParamAry_ref = ();
	my $maxScore = -99999;
	foreach my $minTransfragLen (sort {$a <=> $b} keys %{$optDataForPrintHsh_ref}) {
		foreach my $minTransfragPosCov (sort {$a <=> $b} keys %{$optDataForPrintHsh_ref->{$minTransfragLen}}) {
			foreach my $minTransfragCovPerNt (sort {$a <=> $b} keys %{$optDataForPrintHsh_ref->{$minTransfragLen}{$minTransfragPosCov}}) {
				foreach my $maxTransfragGapTorlerance (sort {$a <=> $b} keys %{$optDataForPrintHsh_ref->{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}}) {
					foreach my $splitWinSize (sort {$a <=> $b} keys %{$optDataForPrintHsh_ref->{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}}) {
						foreach my $splitFold (sort {$a <=> $b} keys %{$optDataForPrintHsh_ref->{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}}) {
							my $score = 0;
							foreach my $dataType (keys %{$optDataForPrintHsh_ref->{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold}}) {
								my $value = $optDataForPrintHsh_ref->{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold}{$dataType};
								$score += $weightFactorHsh_ref->{$dataType}*($value - $dataTypeStatHsh_ref->{$dataType}{'mean'})/$dataTypeStatHsh_ref->{$dataType}{'mean'};
							}
							$optDataForPrintHsh_ref->{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold}{'score'} = $score;
							$optDataForPlotHsh_ref->{'score'}{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold} = $score;
							if ($score > $maxScore) {
								$optimalParamAry_ref = [$minTransfragLen, $minTransfragPosCov, $minTransfragCovPerNt, $maxTransfragGapTorlerance, $splitWinSize, $splitFold];
								$maxScore = $score;
							}
						}
					}
				}
			}
		}
	}
	
	my ($opt_minTransfragLen, $opt_minTransfragPosCov, $opt_minTransfragCovPerNt, $opt_maxTransfragGapTorlerance, $opt_splitWinSize, $opt_splitFold) = @{$optimalParamAry_ref};
	
	my $tmpOptValueHsh_ref = {
		'opt_minTransfragLen' => $opt_minTransfragLen,
		'opt_minTransfragPosCov' => $opt_minTransfragPosCov,
		'opt_minTransfragCovPerNt' => $opt_minTransfragCovPerNt,
		'opt_maxTransfragGapTorlerance' => $opt_maxTransfragGapTorlerance,
		'opt_splitWinSize' => $opt_splitWinSize,
		'opt_splitFold' => $opt_splitFold,
	};
	
	{
		my $LOGFH;
		open $LOGFH, ">", "$resultLogDir/optimized.param.log.txt";
		&printBothFHAndStdout("optimized parameter values:", 10, $LOGFH);#->1495
		&printBothFHAndStdout("$_ = $tmpOptValueHsh_ref->{$_}", 10, $LOGFH) foreach (keys %{$tmpOptValueHsh_ref});#->1495
		&printBothFHAndStdout("optimized score componets:", 10, $LOGFH);#->1495
		foreach my $dataType (sort keys %{$optDataForPrintHsh_ref->{$opt_minTransfragLen}{$opt_minTransfragPosCov}{$opt_minTransfragCovPerNt}{$opt_maxTransfragGapTorlerance}{$opt_splitWinSize}{$opt_splitFold}}) {
			&printBothFHAndStdout("$dataType = $optDataForPrintHsh_ref->{$opt_minTransfragLen}{$opt_minTransfragPosCov}{$opt_minTransfragCovPerNt}{$opt_maxTransfragGapTorlerance}{$opt_splitWinSize}{$opt_splitFold}{$dataType}", 10, $LOGFH);#->1495
		}
		close $LOGFH;
	}
	
	return ($optDataForPrintHsh_ref, $optDataForPlotHsh_ref, $optimalParamAry_ref, $dataTypeStatHsh_ref);
}
sub checkGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: reportStatus|1752
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|163
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: none
#	toCall: &checkGeneInfo($geneInfoHsh_ref);
#	calledInLine: 173
#....................................................................................................................................................#

	my ($geneInfoHsh_ref) = @_;

	&reportStatus("Checking gene categories", 20, "\n");#->1752
	my $ctrgyCountHsh_ref = {};
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		$ctrgyCountHsh_ref->{$ctgry}++;
	}

	foreach my $ctgry (sort keys %{$ctrgyCountHsh_ref}) {
		&reportStatus("Item in $ctgry = $ctrgyCountHsh_ref->{$ctgry}", 20, "\n");#->1752
	}
}
sub checkOverlapAndProximity {
#....................................................................................................................................................#
#	subroutineCategory: range
#	dependOnSub: currentTime|818, generateThreadHshWithRandomCntg|888, reportStatus|1752
#	appearInSub: checkmRNAProximity|633, checkmRNATrnsfrgOverlap|664
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_processInputData|163
#	input: $checkPrxmty, $maxThread, $qryInfoHsh_ref, $qryRngType, $refInfoHsh_ref, $refRngType, $reportExactMatch
#	output: $hitAndPrxmtyByQryHsh_ref, $hitAndPrxmtyByRefHsh_ref
#	toCall: my ($hitAndPrxmtyByRefHsh_ref, $hitAndPrxmtyByQryHsh_ref) = &checkOverlapAndProximity($refInfoHsh_ref, $qryInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType);
#	calledInLine: 657, 682
#....................................................................................................................................................#
	#---incoming variables
	my ($refInfoHsh_ref, $qryInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType) = @_;

	#---outgoing variables
	my $refGeneNumTotal = 0;

	#---make a tmpHsh to contain all cntgs the have either ref and qry
	my $tmpCntgHsh_ref = {};
	my $refCntgHsh_ref = {};
	my $qryCntgHsh_ref = {};

	foreach my $geneID (keys %{$refInfoHsh_ref}) {
		my $cntg = $refInfoHsh_ref->{$geneID}{'cntg'};
		$refGeneNumTotal++;
		$tmpCntgHsh_ref->{$cntg}++;
		$refCntgHsh_ref->{$cntg}{$geneID}++;
	}

	foreach my $geneID (keys %{$qryInfoHsh_ref}) {
		my $cntg = $qryInfoHsh_ref->{$geneID}{'cntg'};
		$tmpCntgHsh_ref->{$cntg}++;
		$qryCntgHsh_ref->{$cntg}{$geneID}++;
	}

	my $totalCntgNum = keys %{$tmpCntgHsh_ref};
	my @cntgAry = keys %{$tmpCntgHsh_ref};

	my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->888
	my $refGeneNumProc :shared = 0;
	my %threadHsh = ();

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\r");#->1752

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781

			sub {
				my ($cntgAry_ref) = @_;

				my $hitAndPrxmtyByRefHsh_InThr_ref = {};
				my $hitAndPrxmtyByQryHsh_InThr_ref = {};

				foreach my $cntg (@{$cntgAry_ref}) {

					#---update the on-screen progress
					#&reportStatus("Finding overlaping on $cntg", 20, "\r");#->1752

					my $tmpPrxmtyByRefHsh_ref = {};
					my $tmpPrxmtyByQryHsh_ref = {};

					if ((exists $qryCntgHsh_ref->{$cntg}) and (exists $refCntgHsh_ref->{$cntg})) {#---if there are both ref and qry can both be found on cntg
						foreach my $refGeneID (keys %{$refCntgHsh_ref->{$cntg}}) {#--- all ftur on the $strnd of $cntg of refGff
							next if not $refInfoHsh_ref->{$refGeneID}{$refRngType};
							my ($refStart, $refEnd) = @{$refInfoHsh_ref->{$refGeneID}{$refRngType}};

							$refGeneNumProc++;

							&reportStatus("$refGeneNumProc of $refGeneNumTotal reference genes checked", 20, "\r");#->1752

							foreach my $qryGeneID (keys %{$qryCntgHsh_ref->{$cntg}}) {#--- all ftur on the $strnd of $cntg of QryGtf
								next if not $qryInfoHsh_ref->{$qryGeneID}{$qryRngType};

								my $samestrnd = "no";
								$samestrnd = "yes" if ($refInfoHsh_ref->{$refGeneID}{'strnd'} eq $qryInfoHsh_ref->{$qryGeneID}{'strnd'});
								my ($qryStart, $qryEnd) = @{$qryInfoHsh_ref->{$qryGeneID}{$qryRngType}};

								my $scene;
								my $ovrlpSize;

								if (($refStart == $qryStart) && ($refEnd == $qryEnd)) {#---scene 0
									$scene = 'exactMatch';
									$ovrlpSize = $qryEnd - $qryStart;

								} elsif (($refStart<=$qryStart)&&($refEnd>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 1
									$scene = 'overlapTail';
									$ovrlpSize = $refEnd - $qryStart;

								} elsif (($refStart>=$qryStart)&&($refStart<=$qryEnd)&&($refEnd>=$qryEnd)) {#---scene 2
									$scene = 'overlapHead';
									$ovrlpSize = $qryEnd - $refStart;

								} elsif (($refStart<=$qryStart)&&($refEnd>=$qryEnd)) {#---scene 3
									$scene = 'cover';
									$ovrlpSize = $qryEnd - $qryStart;

								} elsif (($refStart>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 4
									$scene = 'within';
									$ovrlpSize = $refEnd - $refStart;

								#------Proximity with ref's tail proximal to qry's head
								} elsif (($refEnd<=$qryStart)&&($refEnd<$qryEnd)) {#---scene 5 ---> ref Tail, qry Head

									$scene = 'prxmtyTail';

									if ($checkPrxmty eq "yes") {
										my $tmpPrmxty = $qryStart - $refEnd;
										$tmpPrxmtyByRefHsh_ref->{'XS'}{$refGeneID}{"T"}{$qryGeneID} = $tmpPrmxty;
										$tmpPrxmtyByQryHsh_ref->{'XS'}{$qryGeneID}{"H"}{$refGeneID} = $tmpPrmxty;

										if ($samestrnd eq "yes") {
											$tmpPrxmtyByRefHsh_ref->{'SS'}{$refGeneID}{"T"}{$qryGeneID} = $tmpPrmxty;
											$tmpPrxmtyByQryHsh_ref->{'SS'}{$qryGeneID}{"H"}{$refGeneID} = $tmpPrmxty;
										}
									}

								#------Proximity with ref's head proximal to qry's tail
								} elsif (($refStart>=$qryEnd)&&($refStart>$qryStart)) {#---scene 6 ---> ref Head, qry Tail

									$scene = 'prxmtyHead';

									if ($checkPrxmty eq "yes") {
										my $tmpPrmxty = $refStart - $qryEnd;
										$tmpPrxmtyByRefHsh_ref->{'XS'}{$refGeneID}{"H"}{$qryGeneID} = $tmpPrmxty;
										$tmpPrxmtyByQryHsh_ref->{'XS'}{$qryGeneID}{"T"}{$refGeneID} = $tmpPrmxty;

										if ($samestrnd eq "yes") {
											$tmpPrxmtyByRefHsh_ref->{'SS'}{$refGeneID}{"H"}{$qryGeneID} = $tmpPrmxty;
											$tmpPrxmtyByQryHsh_ref->{'SS'}{$qryGeneID}{"T"}{$refGeneID} = $tmpPrmxty;
										}
									}

								} else {#---BUG! possibly other scene?
									#print "[".&currentTime()."] refStart=$refStart; refEnd=$refEnd; qryStart=$qryStart; qryEnd=$qryEnd\n";#->818
									die "Unexpected overlapping scene between $refGeneID and $qryGeneID. It's a Bug. Program qutting.\n";
								}

								if ($scene ne 'prxmtyTail' and $scene ne 'prxmtyHead' and not ($reportExactMatch eq 'no' and $scene eq 'exactMatch')) {

									@{$hitAndPrxmtyByRefHsh_InThr_ref->{'XS'}{'hit'}{$refGeneID}{$qryGeneID}} = ($scene, $ovrlpSize);
									@{$hitAndPrxmtyByQryHsh_InThr_ref->{'XS'}{'hit'}{$qryGeneID}{$refGeneID}} = ($scene, $ovrlpSize);

									if ($samestrnd eq "yes") {
										@{$hitAndPrxmtyByRefHsh_InThr_ref->{'SS'}{'hit'}{$refGeneID}{$qryGeneID}} = ($scene, $ovrlpSize);
										@{$hitAndPrxmtyByQryHsh_InThr_ref->{'SS'}{'hit'}{$qryGeneID}{$refGeneID}} = ($scene, $ovrlpSize);
									}
								}
							}
						}
					}

					#---find the closest proximity for all refs
					if ($checkPrxmty eq "yes") {
						my $refQryRefHsh_ref = {};

						$refQryRefHsh_ref->{'ref'}{'tmpPrxmtyHsh_ref'} = $tmpPrxmtyByRefHsh_ref;
						$refQryRefHsh_ref->{'ref'}{'cntgHsh_ref'} = $refCntgHsh_ref;
						$refQryRefHsh_ref->{'ref'}{'hitAndPrxmtyHsh_ref'} = $hitAndPrxmtyByRefHsh_InThr_ref;

						$refQryRefHsh_ref->{'qry'}{'tmpPrxmtyHsh_ref'} = $tmpPrxmtyByQryHsh_ref;
						$refQryRefHsh_ref->{'qry'}{'cntgHsh_ref'} = $qryCntgHsh_ref;
						$refQryRefHsh_ref->{'qry'}{'hitAndPrxmtyHsh_ref'} = $hitAndPrxmtyByQryHsh_InThr_ref;

						foreach my $refOrQry ('ref', 'qry') {

							my $cntgHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'cntgHsh_ref'};
							my $tmpPrxmtyHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'tmpPrxmtyHsh_ref'};
							my $hitAndPrxmtyHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'hitAndPrxmtyHsh_ref'};

							foreach my $ftur (keys %{$cntgHsh_ref->{$cntg}}) {
								foreach my $XSOrSS ('XS', 'SS') {
									foreach my $HOrT ('H', 'T') {
										$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{"edge"} = -999 if (not exists $tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT});
										foreach my $otherFtur (sort {$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$a} <=> $tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$b}} keys %{$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}}) {
											@{$hitAndPrxmtyHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} = ($tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$otherFtur}, $otherFtur);
											last; #---sample the smallest only
										}
									}
								}
							}
						}
					}
				}

				return ($hitAndPrxmtyByRefHsh_InThr_ref, $hitAndPrxmtyByQryHsh_InThr_ref);
			}
			,($cntgAry_ref)
		);
	}


	my %tmpTransferThrDataHsh = ();
	$tmpTransferThrDataHsh{'ref'}{'all'} = {};
	$tmpTransferThrDataHsh{'qry'}{'all'} = {};

	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				($tmpTransferThrDataHsh{'ref'}{'thr'}, $tmpTransferThrDataHsh{'qry'}{'thr'}) = $thr->join;
				foreach my $refOrQry (keys %tmpTransferThrDataHsh) {
					my ($allHsh_ref, $thrHsh_ref) = ($tmpTransferThrDataHsh{$refOrQry}{'all'}, $tmpTransferThrDataHsh{$refOrQry}{'thr'});
					foreach my $XSOrSS ('XS', 'SS') {
						foreach my $ftur (keys %{$thrHsh_ref->{$XSOrSS}{'hit'}}) {
							foreach my $hitftur (keys %{$thrHsh_ref->{$XSOrSS}{'hit'}{$ftur}}) {
								@{$allHsh_ref->{$XSOrSS}{'hit'}{$ftur}{$hitftur}} = @{$thrHsh_ref->{$XSOrSS}{'hit'}{$ftur}{$hitftur}};
							}
						}

						if ($thrHsh_ref->{$XSOrSS}{'prxmty'}) {
							foreach my $ftur (keys %{$thrHsh_ref->{$XSOrSS}{'prxmty'}}) {
								foreach my $HOrT (keys %{$thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}}) {
									@{$allHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} = @{$thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} if $thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT};
								}
							}
						}
					}
				}
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	print "\n";

	my $hitAndPrxmtyByRefHsh_ref = $tmpTransferThrDataHsh{'ref'}{'all'};
	my $hitAndPrxmtyByQryHsh_ref = $tmpTransferThrDataHsh{'qry'}{'all'};

	return ($hitAndPrxmtyByRefHsh_ref, $hitAndPrxmtyByQryHsh_ref);
}
sub checkmRNAProximity {
#....................................................................................................................................................#
#	subroutineCategory: range
#	dependOnSub: checkOverlapAndProximity|398, reportStatus|1752
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|163
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref, $mRNAInfoHsh_ref, $maxThread, $resultStorableDir
#	output: $hitAndPrxmtyBymRNAHsh_ref
#	toCall: my ($hitAndPrxmtyBymRNAHsh_ref) = &checkmRNAProximity($resultStorableDir, $geneInfoHsh_ref, $mRNAInfoHsh_ref, $maxThread);
#	calledInLine: 178
#....................................................................................................................................................#
	my ($resultStorableDir, $geneInfoHsh_ref, $mRNAInfoHsh_ref, $maxThread) = @_;

	my $hitAndPrxmtyBymRNAHsh_ref = {};
	my $hitAndPrxmtyBymRNAHshPlsPath = "$resultStorableDir/hitAndPrxmtyBymRNAHsh.pls";
	if (-s $hitAndPrxmtyBymRNAHshPlsPath) {
		&reportStatus("Retrieving hitAndPrxmtyBymRNAHsh", 0,"\n");#->1752
		$hitAndPrxmtyBymRNAHsh_ref = retrieve($hitAndPrxmtyBymRNAHshPlsPath);
	} else {
		&reportStatus("Finding mRNA proximity", 0,"\n");#->1752
		my $checkPrxmty = 'yes';
		my $reportExactMatch = 'no';
		my $refRngType = 'geneRng';
		my $qryRngType = 'geneRng';
		($hitAndPrxmtyBymRNAHsh_ref, undef) = &checkOverlapAndProximity($mRNAInfoHsh_ref, $geneInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType);#->398
		store($hitAndPrxmtyBymRNAHsh_ref, $hitAndPrxmtyBymRNAHshPlsPath);
	}

	return ($hitAndPrxmtyBymRNAHsh_ref);
}
sub checkmRNATrnsfrgOverlap {
#....................................................................................................................................................#
#	subroutineCategory: range
#	dependOnSub: checkOverlapAndProximity|398, reportStatus|1752
#	appearInSub: optimizeTransfrgIdentificationParameters|1254
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_optimizeParameters|185
#	input: $mRNAInfoSubSetHsh_ref, $minOvrlpSize, $trnsfrgInfoHsh_ref
#	output: $antisenseTrnsfrgInfoHsh_ref, $hitToTrnsfrgBymRNAHsh_ref, $hitTomRNAByTrnsfrgHsh_ref, $mRNATrnsfrgOvrlpInfoHsh_ref, $senseTrnsfrgInfoHsh_ref
#	toCall: my ($hitToTrnsfrgBymRNAHsh_ref, $hitTomRNAByTrnsfrgHsh_ref, $mRNATrnsfrgOvrlpInfoHsh_ref, $antisenseTrnsfrgInfoHsh_ref, $senseTrnsfrgInfoHsh_ref) = &checkmRNATrnsfrgOverlap($mRNAInfoSubSetHsh_ref, $trnsfrgInfoHsh_ref, $minOvrlpSize);
#	calledInLine: 1294
#....................................................................................................................................................#

	my ($mRNAInfoSubSetHsh_ref, $trnsfrgInfoHsh_ref, $minOvrlpSize) = @_;

	&reportStatus("Checking overlapping and proximty between mRNA and trnsfrgs", 20, "\n");#->1752

	my $checkPrxmty = 'no';
	my $reportExactMatch = 'yes';
	my ($hitToTrnsfrgBymRNAHsh_ref, $hitTomRNAByTrnsfrgHsh_ref) = &checkOverlapAndProximity($mRNAInfoSubSetHsh_ref, $trnsfrgInfoHsh_ref, $checkPrxmty, $reportExactMatch, 3, 'geneRng', 'geneRng');#->398

	&reportStatus("Pairing mRNA and trnsfrgs", 20, "\n");#->1752
	my $nummRNA = keys %{$mRNAInfoSubSetHsh_ref};
	my $mRNATrnsfrgOvrlpInfoHsh_ref = {};
	my $antisenseTrnsfrgInfoHsh_ref = {};
	my $senseTrnsfrgInfoHsh_ref = {};
	my @paramAry = qw /ovrlpSize covPerNt trnsfrgLen/;

	foreach my $mRNAID (sort keys %{$mRNAInfoSubSetHsh_ref}) {

		#---get gene length
		my ($mRNAStart, $mRNAEnd) = @{$mRNAInfoSubSetHsh_ref->{$mRNAID}{'geneRng'}};
		my $mRNALength = $mRNAEnd - $mRNAStart;
		my $tmpInfoHsh_ref = {};
		my $neighborOvrlpHsh_ref = {};

		#---trnsfrg exists
		if (exists $hitToTrnsfrgBymRNAHsh_ref->{'XS'}{'hit'}{$mRNAID}) {
			foreach my $trnsfrgID (keys %{$hitToTrnsfrgBymRNAHsh_ref->{'XS'}{'hit'}{$mRNAID}}) {
				my $tmpParamValueHsh_ref = {};
				(undef, $tmpParamValueHsh_ref->{'ovrlpSize'}) = @{$hitToTrnsfrgBymRNAHsh_ref->{'XS'}{'hit'}{$mRNAID}{$trnsfrgID}};

				#---skip the trnsfrg if ovrlpSize < $minOvrlpSize;
				next if $tmpParamValueHsh_ref->{'ovrlpSize'} < $minOvrlpSize;

				$tmpParamValueHsh_ref->{'covPerNt'} = $trnsfrgInfoHsh_ref->{$trnsfrgID}{'covPerNt'};
				$tmpParamValueHsh_ref->{'trnsfrgLen'} = $trnsfrgInfoHsh_ref->{$trnsfrgID}{'trnsfrgLen'};

				my $dirtn;
				if ($mRNAInfoSubSetHsh_ref->{$mRNAID}{'strnd'} ne $trnsfrgInfoHsh_ref->{$trnsfrgID}{'strnd'}) {
					$dirtn = 'a';
					%{$antisenseTrnsfrgInfoHsh_ref->{$trnsfrgID}} = %{$trnsfrgInfoHsh_ref->{$trnsfrgID}};
					$antisenseTrnsfrgInfoHsh_ref->{$trnsfrgID}{'neighborOvrlp'} = 'no';
				} else {
					$dirtn = 's';
					%{$senseTrnsfrgInfoHsh_ref->{$trnsfrgID}} = %{$trnsfrgInfoHsh_ref->{$trnsfrgID}};
					$senseTrnsfrgInfoHsh_ref->{$trnsfrgID}{'neighborOvrlp'} = 'no';
				}

				#----push the trsnfrgInfo
				foreach my $param (@paramAry) {
					push @{$tmpInfoHsh_ref->{$dirtn}{$param}}, $tmpParamValueHsh_ref->{$param};
				}

				#----check overlapping
				foreach my $othermRNAID (keys %{$hitTomRNAByTrnsfrgHsh_ref->{'XS'}{'hit'}{$trnsfrgID}}) {
					if ($othermRNAID ne $mRNAID) {
						$neighborOvrlpHsh_ref->{$dirtn}++;
						$antisenseTrnsfrgInfoHsh_ref->{$trnsfrgID}{'neighborOvrlp'} = 'yes' if $dirtn eq 'a';
					}
				}
			}
		}

		#---summary trnsfrg
		foreach my $dirtn ('a', 's') {

			#---zero the ovrlpPct
			$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{'ovrlpPct'} = 0;
			$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{'gapPct'} = 0;
			$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{'trnsfrgNum'} = 0;

			#---record the neighborOvrlp
			$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{'neighborOvrlp'} = 0;
			$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{'neighborOvrlp'} = $neighborOvrlpHsh_ref->{$dirtn} if $neighborOvrlpHsh_ref->{$dirtn};

			#---averge the transfrg param
			foreach my $param (@paramAry) {
				$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{$param} = 0;
				if ($tmpInfoHsh_ref->{$dirtn}{$param}) {
					#---calculate the average for all param
					$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{$param} = sprintf "%.2f", sum(@{$tmpInfoHsh_ref->{$dirtn}{$param}})/@{$tmpInfoHsh_ref->{$dirtn}{$param}};
					#---calculate the ovrlpPct with mRNA
					if ($param eq 'ovrlpSize') {
						$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{'ovrlpPct'} = sprintf "%.2f", 100*sum(@{$tmpInfoHsh_ref->{$dirtn}{'ovrlpSize'}})/$mRNALength;
						$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{'gapPct'} = 100-$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{'ovrlpPct'};
						$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{'trnsfrgNum'} = @{$tmpInfoHsh_ref->{$dirtn}{'ovrlpSize'}};
						$mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{$dirtn}{'avgTrnfrgLen'} = sprintf "%.2f", sum(@{$tmpInfoHsh_ref->{$dirtn}{'trnsfrgLen'}})/@{$tmpInfoHsh_ref->{$dirtn}{'trnsfrgLen'}};
					}
				}
			}
		}
	}

	return ($hitToTrnsfrgBymRNAHsh_ref, $hitTomRNAByTrnsfrgHsh_ref, $mRNATrnsfrgOvrlpInfoHsh_ref, $antisenseTrnsfrgInfoHsh_ref, $senseTrnsfrgInfoHsh_ref);
}
sub classifymRNABasedOnProximity {
#....................................................................................................................................................#
#	subroutineCategory: range, specific
#	dependOnSub: reportStatus|1752
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|163
#	secondaryAppearInSection: >none
#	input: $hitAndPrxmtyBymRNAHsh_ref, $mRNAInfoHsh_ref, $maxPrxmty, $minPrxmty
#	output: $mRNAPrxmtyClassListHsh_ref
#	toCall: my ($mRNAPrxmtyClassListHsh_ref) = &classifymRNABasedOnProximity($hitAndPrxmtyBymRNAHsh_ref, $mRNAInfoHsh_ref, $minPrxmty, $maxPrxmty);
#	calledInLine: 179
#....................................................................................................................................................#
	my ($hitAndPrxmtyBymRNAHsh_ref, $mRNAInfoHsh_ref, $minPrxmty, $maxPrxmty) = @_;

	#---set default valuea
	$minPrxmty = 50 if not $minPrxmty;
	$maxPrxmty = 500 if not $maxPrxmty;

	&reportStatus("Classifying mRNAs based on proximity mRNA", 20,"\n");#->1752

	my %tmpHeadTailHsh = ('H'=>'head', 'T'=>'tail');

	my $mRNAPrxmtyClassListHsh_ref = {};
	foreach my $mRNAID (keys %{$mRNAInfoHsh_ref}) {

		#--exclude any gene that is overlapping with anything
		next if ($hitAndPrxmtyBymRNAHsh_ref->{'XS'}{'hit'}{$mRNAID});

		foreach my $HOrT (keys %{$hitAndPrxmtyBymRNAHsh_ref->{'XS'}{'prxmty'}{$mRNAID}}) {
			my ($prxmtyDist, $prxmtyID) = @{$hitAndPrxmtyBymRNAHsh_ref->{'XS'}{'prxmty'}{$mRNAID}{$HOrT}};
			next if not $mRNAInfoHsh_ref->{$prxmtyID}; #---skip if the proximal feature is not a mRNA
			my $dirtn = 's';
			$dirtn = 'a' if $mRNAInfoHsh_ref->{$mRNAID}{'strnd'} ne $mRNAInfoHsh_ref->{$prxmtyID}{'strnd'};
			if ($prxmtyDist >= $minPrxmty and $prxmtyDist <= $maxPrxmty) {
				$mRNAPrxmtyClassListHsh_ref->{$tmpHeadTailHsh{$HOrT}}{$dirtn}{$mRNAID}++;
			}
		}
	}

	foreach my $headOrTail (keys %{$mRNAPrxmtyClassListHsh_ref}) {
		foreach my $dirtn (keys %{$mRNAPrxmtyClassListHsh_ref->{$headOrTail}}) {
			my $numGene = keys %{$mRNAPrxmtyClassListHsh_ref->{$headOrTail}{$dirtn}};
			&reportStatus("prxmty $headOrTail in dirtn $dirtn = $numGene", 20,"\n");#->1752
		}
	}

	return ($mRNAPrxmtyClassListHsh_ref);
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: >none
#	appearInSub: checkOverlapAndProximity|398, getCtgryGeneInfo|915, getIndivCntgCovPlsPath|949, printBothFHAndStdout|1495, printCMDLogOrFinishMessage|1518, readGFF_oneRNAPerGene|1620, reportStatus|1752
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|65, 5_processInputData|163, 9_finishingTasks|218
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 534, 933, 944, 978, 1512, 1513, 1538, 1541, 1546, 1640, 1768
#....................................................................................................................................................#

	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;

	return $runTime;

}
sub evaluatemRNATrnsfrgOvrlp {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: optimizeTransfrgIdentificationParameters|1254
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_optimizeParameters|185
#	input: $mRNAInfoSubSetHsh_ref, $mRNATrnsfrgOvrlpInfoHsh_ref
#	output: $mRNATrnsfrgOvrlpCountHsh_ref, $mRNAWithSSTrnsfrg
#	toCall: my ($mRNATrnsfrgOvrlpCountHsh_ref, $mRNAWithSSTrnsfrg) = &evaluatemRNATrnsfrgOvrlp($mRNAInfoSubSetHsh_ref, $mRNATrnsfrgOvrlpInfoHsh_ref);
#	calledInLine: 1295
#....................................................................................................................................................#
	my ($mRNAInfoSubSetHsh_ref, $mRNATrnsfrgOvrlpInfoHsh_ref) = @_;

	my $mRNAWithSSTrnsfrg = my $mRNAWithASTrnsfrg = 0;
	my $mRNATrnsfrgOvrlpCountHsh_ref = {};
	$mRNATrnsfrgOvrlpCountHsh_ref->{'sum'}{$_}= 0 foreach qw /NATTrnsfrgNum_AS neighborOvrlp_SS neighborOvrlp_AS mRNATrnsfrgNum_SS mRNACoverPct_SS avgTrnfrgLen_SS avgTrnfrgLen_AS/;
	
	my $totalmRNANum = keys %{$mRNAInfoSubSetHsh_ref};
	foreach my $mRNAID (sort keys %{$mRNAInfoSubSetHsh_ref}) {
		if ($mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{'a'}{'trnsfrgNum'} > 0 and $mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{'a'}{'neighborOvrlp'} == 0) {
			$mRNAWithASTrnsfrg++;
			$mRNATrnsfrgOvrlpCountHsh_ref->{'sum'}{'NATTrnsfrgNum_AS'} += $mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{'a'}{'trnsfrgNum'};
			$mRNATrnsfrgOvrlpCountHsh_ref->{'sum'}{'avgTrnfrgLen_AS'} += $mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{'a'}{'avgTrnfrgLen'};
		}

		if ($mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{'s'}{'trnsfrgNum'} > 0) {
			$mRNAWithSSTrnsfrg++;
			$mRNATrnsfrgOvrlpCountHsh_ref->{'sum'}{'neighborOvrlp_SS'} += $mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{'s'}{'neighborOvrlp'};
			$mRNATrnsfrgOvrlpCountHsh_ref->{'sum'}{'neighborOvrlp_AS'} += $mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{'a'}{'neighborOvrlp'};
			$mRNATrnsfrgOvrlpCountHsh_ref->{'sum'}{'mRNATrnsfrgNum_SS'} += $mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{'s'}{'trnsfrgNum'};
			$mRNATrnsfrgOvrlpCountHsh_ref->{'sum'}{'mRNACoverPct_SS'} += $mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{'s'}{'ovrlpPct'};
			$mRNATrnsfrgOvrlpCountHsh_ref->{'sum'}{'avgTrnfrgLen_SS'} += $mRNATrnsfrgOvrlpInfoHsh_ref->{$mRNAID}{'s'}{'avgTrnfrgLen'};
		}
	}
	
	#---sense 
	foreach (qw /neighborOvrlp_SS neighborOvrlp_AS mRNATrnsfrgNum_SS mRNACoverPct_SS avgTrnfrgLen_SS/) {
		$mRNATrnsfrgOvrlpCountHsh_ref->{'dataToPlot'}{$_} = sprintf "%.5f", $mRNATrnsfrgOvrlpCountHsh_ref->{'sum'}{$_}/$mRNAWithSSTrnsfrg;
	}

	#---anitsense 
	foreach (qw /NATTrnsfrgNum_AS avgTrnfrgLen_AS/) {
		$mRNATrnsfrgOvrlpCountHsh_ref->{'dataToPlot'}{$_} = sprintf "%.5f", $mRNATrnsfrgOvrlpCountHsh_ref->{'sum'}{$_}/$mRNAWithASTrnsfrg;
	}
	
	$mRNATrnsfrgOvrlpCountHsh_ref->{'dataToPlot'}{'pctWithTrnsfrg_AS'} = sprintf "%.5f", 100*$mRNAWithASTrnsfrg/$totalmRNANum;
	$mRNATrnsfrgOvrlpCountHsh_ref->{'dataToPlot'}{'pctWithTrnsfrg_SS'} = sprintf "%.5f", 100*$mRNAWithSSTrnsfrg/$totalmRNANum;

	return ($mRNATrnsfrgOvrlpCountHsh_ref, $mRNAWithSSTrnsfrg);
}
sub generateThreadHshWithRandomCntg {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: checkOverlapAndProximity|398, identifyTransfrag|1078
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $cntgAry_ref, $threadToSpawn
#	output: $randCntgInThreadHsh_ref
#	toCall: my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($threadToSpawn, $cntgAry_ref);
#	calledInLine: 436, 1097
#....................................................................................................................................................#

	my ($threadToSpawn, $cntgAry_ref) = @_;

	my @shuffleCntgAry = shuffle(@{$cntgAry_ref});
	my $threadNum = 1;
	my $randCntgInThreadHsh_ref = {};
	foreach my $cntg (@{$cntgAry_ref}) {
		$threadNum = 1 if $threadNum > $threadToSpawn;
		push @{$randCntgInThreadHsh_ref->{$threadNum}}, $cntg;
		$threadNum++;
	}

	return $randCntgInThreadHsh_ref;

}
sub getCtgryGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|818
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|163
#	secondaryAppearInSection: >none
#	input: $ctgryAry_ref, $geneInfoHsh_ref
#	output: $geneCtgryInfoHsh_ref
#	toCall: my ($geneCtgryInfoHsh_ref) = &getCtgryGeneInfo($geneInfoHsh_ref, $ctgryAry_ref);
#	calledInLine: 177
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $ctgryAry_ref) = @_;

	my $geneCtgryInfoHsh_ref = {};

	my $ctgryStr = join ",", @{$ctgryAry_ref};

	print "[".&currentTime()."] Filtering GFF on cgtry $ctgryStr.\n";#->818

	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		if (grep /^$ctgry$/, @{$ctgryAry_ref}) {
			%{$geneCtgryInfoHsh_ref->{$geneID}} = %{$geneInfoHsh_ref->{$geneID}};
		}
	}

	my $numGene = keys %{$geneCtgryInfoHsh_ref};

	print "[".&currentTime()."] $numGene gene filtered on cgtry $ctgryStr.\n";#->818

	return $geneCtgryInfoHsh_ref;
}
sub getIndivCntgCovPlsPath {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: currentTime|818
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|163
#	secondaryAppearInSection: >none
#	input: $cntgCovPlsIndexPath
#	output: $cntgCovInPlsPathHsh_ref
#	toCall: my ($cntgCovInPlsPathHsh_ref) = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);
#	calledInLine: 169, 961
#....................................................................................................................................................#

	#my $cntgCovInPlsPathHsh_ref = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);#->949
	my ($cntgCovPlsIndexPath) = @_;

	$cntgCovPlsIndexPath =~ s/\.gz$//;

	my $cntgCovInPlsPathHsh_ref = {};

	system ("gzip -fd $cntgCovPlsIndexPath.gz") if -s "$cntgCovPlsIndexPath.gz";
	my %plsIndexHsh = %{retrieve($cntgCovPlsIndexPath)};

	my (undef, $cntgCovStroableDir, undef) = fileparse($cntgCovPlsIndexPath, qr/\.[^.]*/);
	foreach my $cntg (keys %plsIndexHsh) {
		my $cntgCovPlsPath = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
		die "cntgCovPlsPath $cntgCovPlsPath is invalid\n" if ((not -s $cntgCovPlsPath) and (not -s $cntgCovPlsPath.".gz"));
		$cntgCovInPlsPathHsh_ref->{$cntg} = $cntgCovPlsPath;
	}
	my $numCntg = keys %{$cntgCovInPlsPathHsh_ref};
	print "[".&currentTime()."] pls path of $numCntg contig stored.\n";#->818

	return $cntgCovInPlsPathHsh_ref;
}
sub getTransfragUsingOptimalParameters {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: identifyTransfrag|1078
#	appearInSub: >none
#	primaryAppearInSection: 7_identifyTransfragsUsingOptParam|197
#	secondaryAppearInSection: >none
#	input: $maxThread, $optimalParamAry_ref, $queryCovPlsPathHsh_ref, $resultStorableDir
#	output: $trnsfrgInfoHsh_ref
#	toCall: my ($trnsfrgInfoHsh_ref) = &getTransfragUsingOptimalParameters($optimalParamAry_ref, $queryCovPlsPathHsh_ref, $maxThread, $resultStorableDir);
#	calledInLine: 202
#....................................................................................................................................................#
	my ($optimalParamAry_ref, $queryCovPlsPathHsh_ref, $maxThread, $resultStorableDir) = @_;
	
	my ($minTransfragLen, $minTransfragPosCov, $minTransfragCovPerNt, $maxTransfragGapTorlerance, $splitWinSize, $splitFold) = @{$optimalParamAry_ref};
	my $keepAllSplit = 'no';
	my $minCovPerNtInWindow = 10;
	my $splitStepSize = 1;
	my ($trnsfrgInfoHsh_ref) = &identifyTransfrag($queryCovPlsPathHsh_ref, $minTransfragLen, $minTransfragPosCov, $minTransfragCovPerNt, $maxTransfragGapTorlerance, $splitWinSize, $splitFold, $keepAllSplit, $minCovPerNtInWindow, $splitStepSize, $maxThread);#->1078

	my $trnsfrgInfoHshPlsPath = "$resultStorableDir/trnsfrgInfoHsh.pls";
	store($trnsfrgInfoHsh_ref, $trnsfrgInfoHshPlsPath);

	return ($trnsfrgInfoHsh_ref);
}
sub getmRNASubset {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|1752
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|163
#	secondaryAppearInSection: >none
#	input: $mRNAInfoHsh_ref, $mRNAPrxmtyClassListHsh_ref, $queryCovPlsPathHsh_ref
#	output: $mRNAInfoSubSetHsh_ref, $queryCovSubsetPlsPathHsh_ref
#	toCall: my ($mRNAInfoSubSetHsh_ref, $queryCovSubsetPlsPathHsh_ref) = &getmRNASubset($mRNAInfoHsh_ref, $mRNAPrxmtyClassListHsh_ref, $queryCovPlsPathHsh_ref);
#	calledInLine: 180
#....................................................................................................................................................#
	my ($mRNAInfoHsh_ref, $mRNAPrxmtyClassListHsh_ref, $queryCovPlsPathHsh_ref) = @_;

	my $mRNAInfoSubSetHsh_ref = {};
	my $queryCovSubsetPlsPathHsh_ref = {};

	#---get only genes with neighbor gene that is in dirtn sense
	foreach my $HOrT (keys %{$mRNAPrxmtyClassListHsh_ref}) {
		foreach my $mRNAID (keys %{$mRNAPrxmtyClassListHsh_ref->{$HOrT}{'s'}}) {
			$mRNAInfoSubSetHsh_ref->{$mRNAID} = $mRNAInfoHsh_ref->{$mRNAID} if not $mRNAInfoSubSetHsh_ref->{$mRNAID};
			$queryCovSubsetPlsPathHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}} = $queryCovPlsPathHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}} if not $queryCovSubsetPlsPathHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}};
		}
	}

	my $numGeneInSubSet = keys %{$mRNAInfoSubSetHsh_ref};
	my $numCntgInSubSet = keys %{$queryCovSubsetPlsPathHsh_ref};
	&reportStatus("$numGeneInSubSet mRNA in $numCntgInSubSet cntg were selected as the subSet for parameter optimization", 10, "\n");#->1752

	return ($mRNAInfoSubSetHsh_ref, $queryCovSubsetPlsPathHsh_ref);
}
sub ggplotXYLineThreeFactors {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: plotTransfrgIdentificationParameters|1424
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_optimizeParameters|185
#	input: $RScriptPath, $YAxis, $dataPath, $extraArg, $extraStatment, $factorA, $factorB, $factorC, $height, $logPath, $pdfPath, $plotDataHsh_ref, $width
#	output: 
#	toCall: &ggplotXYLineThreeFactors($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $factorA, $factorB, $factorC, $YAxis, $extraStatment, $extraArg, $height, $width);
#	calledInLine: 1485
#....................................................................................................................................................#

	my ($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $factorA, $factorB, $factorC, $YAxis, $extraStatment, $extraArg, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA join "", (join "\t", ($factorA, $factorB, $factorC, $YAxis)), "\n";
	foreach my $AVal (sort {$a <=> $b} keys %{$plotDataHsh_ref}) {
		foreach my $BVal (sort {$a <=> $b} keys %{$plotDataHsh_ref->{$AVal}}) {
			foreach my $CVal (sort {$a <=> $b} keys %{$plotDataHsh_ref->{$AVal}{$BVal}}) {
				my $YVal = $plotDataHsh_ref->{$AVal}{$BVal}{$CVal};
				print PLOTDATA join "", (join "\t", ($AVal, $BVal, $CVal, $YVal)), "\n";
			}
		}
	}
	close PLOTDATA;

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "$extraStatment"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$factorA, y=$YAxis)) + geom_point(aes(colour = $factorB), size = 5, alpha = 0.8) + scale_colour_gradient(low = \"blue\", high = \"red\") + facet_grid(~$factorC) $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath");
	
	return ();
}
sub identifyTransfrag {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: generateThreadHshWithRandomCntg|888, reportStatus|1752, scanDramaticCoverageChange|1773
#	appearInSub: getTransfragUsingOptimalParameters|983, optimizeTransfrgIdentificationParameters|1254
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_optimizeParameters|185, 7_identifyTransfragsUsingOptParam|197
#	input: $keepAllSplit, $maxThread, $maxTransfragGapTorlerance, $minCovPerNtInWindow, $minTransfragCovPerNt, $minTransfragLen, $minTransfragPosCov, $queryCovPlsPathHsh_ref, $splitFold, $splitStepSize, $splitWinSize
#	output: $trnsfrgInfoHsh_ref
#	toCall: my ($trnsfrgInfoHsh_ref) = &identifyTransfrag($queryCovPlsPathHsh_ref, $minTransfragLen, $minTransfragPosCov, $minTransfragCovPerNt, $maxTransfragGapTorlerance, $splitWinSize, $splitFold, $keepAllSplit, $minCovPerNtInWindow, $splitStepSize, $maxThread);
#	calledInLine: 1000, 1293
#....................................................................................................................................................#

	my ($queryCovPlsPathHsh_ref, $minTransfragLen, $minTransfragPosCov, $minTransfragCovPerNt, $maxTransfragGapTorlerance, $splitWinSize, $splitFold, $keepAllSplit, $minCovPerNtInWindow, $splitStepSize, $maxThread) = @_;

	my $trnsfrgInfoHsh_ref = {};

	my $totalCntgNum = keys %{$queryCovPlsPathHsh_ref};
	my @cntgAry = keys %{$queryCovPlsPathHsh_ref};

	my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->888
	my $numCntgProc :shared = 0;
	my $numTrnsfrg :shared = 0;
	my %threadHsh = ();

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\r");#->1752

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781

			sub {
				my ($cntgAry_ref) = @_;

				my $trnsfrgInfoHsh_inThr_ref = {};

				foreach my $cntg (@{$cntgAry_ref}) {
					my $cntgCovPlsPath = $queryCovPlsPathHsh_ref->{$cntg};
					$numCntgProc++;

					system ("gzip -df $cntgCovPlsPath.gz") if (-s "$cntgCovPlsPath.gz");
					my $cntgCovAry_ref = retrieve($cntgCovPlsPath);

					#----define tmp variables
					my $tmpGapLenHsh_ref = {};
					my $tmpBreakHsh_ref = {};
					my $tmpTrnsfrgPosHsh_ref = {};
					my $tmpTrnsfrgCovHsh_ref = {};

					foreach my $strnd ('+', '-') {
						$tmpBreakHsh_ref->{$strnd} = 'yes';
						$tmpGapLenHsh_ref->{$strnd} = 0;
					}

					#----go through cntg position by position
					my $numTrnsfrgOnCntg = 0;

					foreach my $i (0..$#{$cntgCovAry_ref}) {
						my $pos = $i + 1;
						my $tmpPosCovHsh_ref = {};
						$tmpPosCovHsh_ref->{'+'} = 0;
						$tmpPosCovHsh_ref->{'-'} = 0;
						($tmpPosCovHsh_ref->{'+'}, $tmpPosCovHsh_ref->{'-'}) = split /,/, $cntgCovAry_ref->[$i] if $cntgCovAry_ref->[$i];

						foreach my $strnd ('+', '-') {#---simultaneously run thr two strnd

							#-------check the status of the position to determine it's a break or not
							if ($tmpPosCovHsh_ref->{$strnd} >= $minTransfragPosCov) {#----over minimum cov
								push @{$tmpTrnsfrgPosHsh_ref->{$strnd}}, $pos;
								push @{$tmpTrnsfrgCovHsh_ref->{$strnd}}, $tmpPosCovHsh_ref->{$strnd};

								#---reset the gap and break
								$tmpGapLenHsh_ref->{$strnd} = 0;
								$tmpBreakHsh_ref->{$strnd} = 'no';

							} elsif ($tmpBreakHsh_ref->{$strnd} eq 'no') {#----within a gap but not inside a long break
								$tmpGapLenHsh_ref->{$strnd}++;
								#----within a gap but torlerated
								if ($tmpGapLenHsh_ref->{$strnd} <= $maxTransfragGapTorlerance) {
									push @{$tmpTrnsfrgPosHsh_ref->{$strnd}}, $pos;
									push @{$tmpTrnsfrgCovHsh_ref->{$strnd}}, $tmpPosCovHsh_ref->{$strnd};

								} else {#----within a gap and exceeded
									$tmpBreakHsh_ref->{$strnd} = 'yes';
								}
							}

							#-----if this is the end of the cntg it will be treated as a long break
							$tmpBreakHsh_ref->{$strnd} = 'yes' if $i == $#{$cntgCovAry_ref};

							#------Inside the break and has something in $tmpTrnsfrgPosHsh_ref
							if ($tmpBreakHsh_ref->{$strnd} eq 'yes' and $tmpTrnsfrgPosHsh_ref->{$strnd}) {
								#-----$trnsfrgLen = number of positions recorded minus the gap len
								my $trnsfrgLen = @{$tmpTrnsfrgPosHsh_ref->{$strnd}} - $tmpGapLenHsh_ref->{$strnd};

								if ($trnsfrgLen >= $minTransfragLen) {

									$numTrnsfrgOnCntg++;

									#----remove the trailing gaps
									for my $j (1..$tmpGapLenHsh_ref->{$strnd}) {
										pop @{$tmpTrnsfrgPosHsh_ref->{$strnd}};
										pop @{$tmpTrnsfrgCovHsh_ref->{$strnd}};
									}

									my $splitTrnsfrgRngHsh_ref = &scanDramaticCoverageChange($tmpTrnsfrgPosHsh_ref, $tmpTrnsfrgCovHsh_ref, $strnd, $splitWinSize, $splitFold, $minCovPerNtInWindow, $splitStepSize);#->1773

									my $splitNum = 0;
									foreach my $startIndex (sort keys %{$splitTrnsfrgRngHsh_ref}) {
										my $endIndex = $splitTrnsfrgRngHsh_ref->{$startIndex};
										my $covSum = sprintf "%.1f", sum(@{$tmpTrnsfrgCovHsh_ref->{$strnd}}[$startIndex..$endIndex]);
										my $trnsfrgLen = $endIndex-$startIndex;
										my $covPerNt = sprintf "%.1f", $covSum/$trnsfrgLen;
										if ($covPerNt >= $minTransfragCovPerNt) {
											next if $trnsfrgLen < $minTransfragLen and $keepAllSplit eq 'no';
											my $direction = "F"; $direction = "R" if $strnd eq '-';
											my $trnsfrgID = join "_", ($cntg, "$numTrnsfrgOnCntg.$splitNum", $direction, $covPerNt, $trnsfrgLen);
											my $description = "covPerNt[$covPerNt]";
											my $trnsfrgStart = ${$tmpTrnsfrgPosHsh_ref->{$strnd}}[$startIndex];
											my $trnsfrgEnd = ${$tmpTrnsfrgPosHsh_ref->{$strnd}}[$endIndex];
											$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'cntg'} = $cntg;
											$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'covPerNt'} = $covPerNt;
											$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'trnsfrgLen'} = $trnsfrgLen;
											$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'trnsfrgStart'} = $trnsfrgStart;
											$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'trnsfrgEnd'} = $trnsfrgEnd;
											$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'strnd'} = $strnd;
											$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'description'} = $description;
											@{$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'geneRng'}} = ($trnsfrgStart, $trnsfrgEnd);
											$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'ctgry'} = $cntg;
											$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'RNAID'} = "rna_".$trnsfrgID;
											@{$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'RNARng'}} = ($trnsfrgStart, $trnsfrgEnd);
											@{$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'geneRng'}} = ($trnsfrgStart, $trnsfrgEnd);
											push @{$trnsfrgInfoHsh_inThr_ref->{$trnsfrgID}{'exonRng'}}, ($trnsfrgStart, $trnsfrgEnd);
											$numTrnsfrg++;
										}
										$splitNum++;
									}

									&reportStatus("$numTrnsfrg trnsfrgs identified in $numCntgProc cntg", 20, "\r");#->1752

								}#---end of if ($trnsfrgLen >= $minTransfragLen)

								$tmpGapLenHsh_ref->{$strnd} = 0; #---reset gap
								delete $tmpTrnsfrgPosHsh_ref->{$strnd}; #---reset pos
								delete $tmpTrnsfrgCovHsh_ref->{$strnd}; #---reset cov

							}#---end of if ($tmpBreakHsh_ref->{$strnd} eq 'yes' and $tmpTrnsfrgPosHsh_ref->{$strnd})
						}#---end of foreach my $strnd ('+', '-')
					}#---end of foreach my $i (0..$#{$cntgCovAry_ref})
				}#---end of foreach my $cntg (sort keys %{$queryCovPlsPathHsh_ref})

				return ($trnsfrgInfoHsh_inThr_ref);
			}#---end of sub
			,($cntgAry_ref)
		);#---end of threads->new
	}#---end of foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref})

	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				my ($trnsfrgInfoHsh_inThr_ref) = $thr->join;
				%{$trnsfrgInfoHsh_ref->{$_}} = %{$trnsfrgInfoHsh_inThr_ref->{$_}} foreach (keys %{$trnsfrgInfoHsh_inThr_ref});
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	my $storedTrnsfrg = keys %{$trnsfrgInfoHsh_ref};
	&reportStatus("$storedTrnsfrg transfrags were stored", 20, "\n");#->1752

	return ($trnsfrgInfoHsh_ref);
}
sub optimizeTransfrgIdentificationParameters {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkmRNATrnsfrgOverlap|664, evaluatemRNATrnsfrgOvrlp|837, identifyTransfrag|1078, reportStatus|1752
#	appearInSub: >none
#	primaryAppearInSection: 6_optimizeParameters|185
#	secondaryAppearInSection: >none
#	input: $mRNAInfoSubSetHsh_ref, $maxThread, $maxTransfragGapTorleranceAry_ref, $minTransfragCovPerNtAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $queryCovSubsetPlsPathHsh_ref, $resultStorableDir, $splitFoldAry_ref, $splitWinSizeAry_ref
#	output: $optTransfrgParamDataHsh_ref
#	toCall: my ($optTransfrgParamDataHsh_ref) = &optimizeTransfrgIdentificationParameters($queryCovSubsetPlsPathHsh_ref, $mRNAInfoSubSetHsh_ref, $maxTransfragGapTorleranceAry_ref, $splitWinSizeAry_ref, $splitFoldAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $minTransfragCovPerNtAry_ref, $maxThread, $resultStorableDir);
#	calledInLine: 190
#....................................................................................................................................................#
	my ($queryCovSubsetPlsPathHsh_ref, $mRNAInfoSubSetHsh_ref, $maxTransfragGapTorleranceAry_ref, $splitWinSizeAry_ref, $splitFoldAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $minTransfragCovPerNtAry_ref, $maxThread, $resultStorableDir) = @_;

	my $optTransfrgParamDataHsh_ref = {};
	my $optTransfrgParamDataHshPlsPath = "$resultStorableDir/optTransfrgParamDataHsh.pls";
	$optTransfrgParamDataHsh_ref = retrieve($optTransfrgParamDataHshPlsPath) if -s $optTransfrgParamDataHshPlsPath;

	my $keepAllSplit = 'yes';
	my $minCovPerNtInWindow = 10;
	my $minOvrlpSize = 5;
	my $combinationProc = 0;
	my $numCombination = @{$maxTransfragGapTorleranceAry_ref}*@{$splitWinSizeAry_ref}*@{$splitFoldAry_ref}*@{$minTransfragLenAry_ref}*@{$minTransfragPosCovAry_ref}*@{$minTransfragCovPerNtAry_ref};
	foreach my $minTransfragLen (sort {$a <=> $b} @{$minTransfragLenAry_ref}) {
		foreach my $minTransfragPosCov (sort {$a <=> $b} @{$minTransfragPosCovAry_ref}) {
			foreach my $minTransfragCovPerNt (sort {$a <=> $b} @{$minTransfragCovPerNtAry_ref}) {
				foreach my $maxTransfragGapTorlerance (sort {$a <=> $b} @{$maxTransfragGapTorleranceAry_ref}) {
					foreach my $splitWinSize (sort {$a <=> $b} @{$splitWinSizeAry_ref}) {
						foreach my $splitFold (sort {$a <=> $b} @{$splitFoldAry_ref}) {
							$combinationProc++;
							my $splitStepSize = $splitWinSize/4;

							my $paramTag = "ML$minTransfragLen.MC$minTransfragPosCov.CN$minTransfragCovPerNt.GP$maxTransfragGapTorlerance.SW$splitWinSize.SF$splitFold";
							&reportStatus("Processing $combinationProc of $numCombination combinations: $paramTag", 10, "\n");#->1752

							if ($optTransfrgParamDataHsh_ref->{'neighborOvrlp_SS'}{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold}) {
								&reportStatus("Results found. Skipping", 10, "\n");#->1752
								next;
							}
							
							my ($trnsfrgInfoHsh_ref) = &identifyTransfrag($queryCovSubsetPlsPathHsh_ref, $minTransfragLen, $minTransfragPosCov, $minTransfragCovPerNt, $maxTransfragGapTorlerance, $splitWinSize, $splitFold, $keepAllSplit, $minCovPerNtInWindow, $splitStepSize, $maxThread);#->1078
							my ($hitToTrnsfrgBymRNAHsh_ref, $hitTomRNAByTrnsfrgHsh_ref, $mRNATrnsfrgOvrlpInfoHsh_ref, $antisenseTrnsfrgInfoHsh_ref, $senseTrnsfrgInfoHsh_ref) = &checkmRNATrnsfrgOverlap($mRNAInfoSubSetHsh_ref, $trnsfrgInfoHsh_ref, $minOvrlpSize);#->664
							my ($mRNATrnsfrgOvrlpCountHsh_ref, $mRNAWithSSTrnsfrg) = &evaluatemRNATrnsfrgOvrlp($mRNAInfoSubSetHsh_ref, $mRNATrnsfrgOvrlpInfoHsh_ref);#->837
							foreach (keys %{$mRNATrnsfrgOvrlpCountHsh_ref->{'dataToPlot'}}) {
								$optTransfrgParamDataHsh_ref->{$_}{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold} = $mRNATrnsfrgOvrlpCountHsh_ref->{'dataToPlot'}{$_};
							}
						}
					}
				}
			}
		}
	}
	
	store($optTransfrgParamDataHsh_ref, $optTransfrgParamDataHshPlsPath);

	return ($optTransfrgParamDataHsh_ref);
}
sub outputFinalTransfragGff {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: printGFF_oneRNAPerGene_chooseStrnd_filter|1551
#	appearInSub: >none
#	primaryAppearInSection: 8_outputFiles|207
#	secondaryAppearInSection: >none
#	input: $IGVGffTrackHsh_ref, $resultGffDir, $trnsfrgInfoHsh_ref
#	output: 
#	toCall: &outputFinalTransfragGff($resultGffDir, $trnsfrgInfoHsh_ref, $IGVGffTrackHsh_ref);
#	calledInLine: 212
#....................................................................................................................................................#
	my ($resultGffDir, $trnsfrgInfoHsh_ref, $IGVGffTrackHsh_ref) = @_;
	
	my $trnsfragGFFPathHsh_ref = {
		'plus' => "$resultGffDir/trnsfrag.plus.gff",
		'minus' => "$resultGffDir/trnsfrag.minus.gff",
	};
	
	foreach my $plusOrMinus (keys %{$trnsfragGFFPathHsh_ref}) {
		my $trnsfragGFFPath = $trnsfragGFFPathHsh_ref->{$plusOrMinus};
		my $geneInfoHsh_ref = $trnsfrgInfoHsh_ref;
		my $outGFFPath = $trnsfragGFFPath;
		my $gffGeneLineOnly = 'yes';
		my $strndInWord = $plusOrMinus;
		my $filter = undef;
		&printGFF_oneRNAPerGene_chooseStrnd_filter($geneInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, $filter);#->1551
		$IGVGffTrackHsh_ref->{'trnsfrg'}{$plusOrMinus} = $trnsfragGFFPath;
	}

	return ();
}
sub outputIGVXML {
#....................................................................................................................................................#
#	subroutineCategory: XML
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 8_outputFiles|207
#	secondaryAppearInSection: >none
#	input: $IGVGenomeID, $IGVGffTrackHsh_ref, $XMLPath, $basicIGVTrackHsh_ref
#	output: none
#	toCall: &outputIGVXML($XMLPath, $basicIGVTrackHsh_ref, $IGVGenomeID, $IGVGffTrackHsh_ref);
#	calledInLine: 213
#....................................................................................................................................................#

	my ($XMLPath, $basicIGVTrackHsh_ref, $IGVGenomeID, $IGVGffTrackHsh_ref) = @_;

	my $colorHsh_ref = {};
	$colorHsh_ref->{'plus'} = '255,153,153';
	$colorHsh_ref->{'minus'} = '153,153,255';

	open (XML, ">", $XMLPath);
	printf XML "%0s", "<\?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"\?>\n";
	printf XML "%0s", "<Session genome=\"$IGVGenomeID\" version=\"5\">\n";
	printf XML "%4s", "<Resources>\n";

	#---Resource for basic tracks
	foreach my $trackName (keys %{$basicIGVTrackHsh_ref}) {
		my $filePath = $basicIGVTrackHsh_ref->{$trackName};
		printf XML "%8s", "<Resource path=\"$filePath\"\/>\n";
	}

	#---Resource for gff tracks
	foreach my $trackName (keys %{$IGVGffTrackHsh_ref}) {
		foreach my $plusOrMinus (keys %{$IGVGffTrackHsh_ref->{$trackName}}) {
			my $filePath = $IGVGffTrackHsh_ref->{$trackName}{$plusOrMinus};
			printf XML "%8s", "<Resource path=\"$filePath\"\/>\n";
		}
	}

	printf XML "%4s", "</Resources>\n";
	printf XML "%4s", "<Panel height=\"900\" name=\"DataPanel\" width=\"1200\">\n";

	printf XML "%8s", "<Track colorScale=\"ContinuousColorScale;20.0;10.0;25.0;35.0;153,153,255;255,255,255;255,153,153\" fontSize=\"10\" id=\"$basicIGVTrackHsh_ref->{'GCPath'}\" name=\"GC\" normalize=\"false\" renderer=\"HEATMAP\" height=\"20\" sortable=\"true\" visible=\"true\" windowFunction=\"none\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track colorScale=\"ContinuousColorScale;2.0;1.0;3.0;15.0;255,255,255;153,153,255;255,153,153\" fontSize=\"10\" id=\"$basicIGVTrackHsh_ref->{'reptvPath'}\" name=\"Repetitiveness\" normalize=\"false\" renderer=\"HEATMAP\" height=\"20\" sortable=\"true\" visible=\"true\" windowFunction=\"none\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" colorScale=\"ContinuousColorScale;0.0;17.0;255,255,255;0,0,178\" displayMode=\"SQUISHED\"  fontSize=\"10\" height=\"20\" id=\"EHI_v13_genes\" name=\"bonaFide mRNA\" renderer=\"BASIC_FEATURE\" sortable=\"false\" visible=\"true\" windowFunction=\"count\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" colorScale=\"ContinuousColorScale;0.0;40.0;255,255,255;0,0,178\" displayMode=\"COLLAPSED\"  fontSize=\"10\" height=\"20\" id=\"$basicIGVTrackHsh_ref->{'gffPath'}\" name=\"All Genomic Features\" renderer=\"GENE_TRACK\" sortable=\"false\" visible=\"true\" windowFunction=\"count\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"true\" color=\"255,153,153\" displayMode=\"COLLAPSED\" fontSize=\"10\" height=\"100\" id=\"$basicIGVTrackHsh_ref->{'plusCovTDFPath'}\" name=\"Corrected Plus Cov\" normalize=\"false\" renderer=\"BAR_CHART\" sortable=\"true\" visible=\"true\" windowFunction=\"mean\">\n";
	printf XML "%12s", "<DataRange type=\"LOG\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"true\" color=\"153,153,255\" displayMode=\"COLLAPSED\" fontSize=\"10\" height=\"100\" id=\"$basicIGVTrackHsh_ref->{'minusCovTDFPath'}\" name=\"Corrected Minus Cov\" normalize=\"false\" renderer=\"BAR_CHART\" sortable=\"true\" visible=\"true\" windowFunction=\"mean\">\n";
	printf XML "%12s", "<DataRange type=\"LOG\"/>\n";
	printf XML "%8s", "</Track>\n";

	foreach my $trackName (keys %{$IGVGffTrackHsh_ref}) {
		foreach my $plusOrMinus (keys %{$IGVGffTrackHsh_ref->{$trackName}}) {
			my $filePath = $IGVGffTrackHsh_ref->{$trackName}{$plusOrMinus};
			printf XML "%8s", "<Track autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"$colorHsh_ref->{$plusOrMinus}\" displayMode=\"EXPANDED\"  fontSize=\"10\" height=\"20\" id=\"$filePath\" name=\"$trackName\_$plusOrMinus\" renderer=\"GENE_TRACK\" sortable=\"false\" visible=\"true\" windowFunction=\"count\">\n";
			printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
			printf XML "%8s", "</Track>\n";
		}
	}

	printf XML "%4s", "</Panel>\n";

	printf XML "%0s", "</Session>\n";

	close XML;

}
sub plotTransfrgIdentificationParameters {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotXYLineThreeFactors|1039, reportStatus|1752
#	appearInSub: >none
#	primaryAppearInSection: 6_optimizeParameters|185
#	secondaryAppearInSection: >none
#	input: $dataTypeStatHsh_ref, $ggplotDirHsh_ref, $maxTransfragGapTorleranceAry_ref, $minTransfragCovPerNtAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $optDataForPlotHsh_ref, $optimalParamAry_ref, $splitFoldAry_ref, $splitWinSizeAry_ref
#	output: 
#	toCall: &plotTransfrgIdentificationParameters($optDataForPlotHsh_ref, $ggplotDirHsh_ref, $optimalParamAry_ref, $maxTransfragGapTorleranceAry_ref, $splitWinSizeAry_ref, $splitFoldAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $minTransfragCovPerNtAry_ref, $dataTypeStatHsh_ref);
#	calledInLine: 192
#....................................................................................................................................................#
	my ($optDataForPlotHsh_ref, $ggplotDirHsh_ref, $optimalParamAry_ref, $maxTransfragGapTorleranceAry_ref, $splitWinSizeAry_ref, $splitFoldAry_ref, $minTransfragLenAry_ref, $minTransfragPosCovAry_ref, $minTransfragCovPerNtAry_ref, $dataTypeStatHsh_ref) = @_;
	
	my ($opt_minTransfragLen, $opt_minTransfragPosCov, $opt_minTransfragCovPerNt, $opt_maxTransfragGapTorlerance, $opt_splitWinSize, $opt_splitFold) = @{$optimalParamAry_ref};
	
	foreach my $dataType (keys %{$optDataForPlotHsh_ref}) {
		my $dirTag = $dataType;
		system ("mkdir -pm 777 $ggplotDirHsh_ref->{$_}/$dirTag/") foreach (keys %{$ggplotDirHsh_ref});
		foreach my $minTransfragLen (sort {$a <=> $b} @{$minTransfragLenAry_ref}) {
			foreach my $minTransfragPosCov (sort {$a <=> $b} @{$minTransfragPosCovAry_ref}) {
				foreach my $minTransfragCovPerNt (sort {$a <=> $b} @{$minTransfragCovPerNtAry_ref}) {

					my $transfragLenCovTag = "ML$minTransfragLen\_MC$minTransfragPosCov\_CN$minTransfragCovPerNt";
					my $plotDataHsh_ref = {};

					&reportStatus("Plotting $dataType for $transfragLenCovTag", 10, "\n");#->1752
					
					foreach my $maxTransfragGapTorlerance (sort {$a <=> $b} @{$maxTransfragGapTorleranceAry_ref}) {
						foreach my $splitWinSize (sort {$a <=> $b} @{$splitWinSizeAry_ref}) {
							foreach my $splitFold (sort {$a <=> $b} @{$splitFoldAry_ref}) {

								my $value = $optDataForPlotHsh_ref->{$dataType}{$minTransfragLen}{$minTransfragPosCov}{$minTransfragCovPerNt}{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold};
								$plotDataHsh_ref->{$maxTransfragGapTorlerance}{$splitWinSize}{$splitFold} = $value;
							}
						}
					}
					
					{
						#---add the optimal value to the plot if the transfragLenCov parameter is optimal
						my $extraStatment = my $extraArg = '';
						if ($opt_minTransfragLen eq $minTransfragLen and $opt_minTransfragPosCov eq $minTransfragPosCov and $opt_minTransfragCovPerNt eq $minTransfragCovPerNt) {
							my $optimalValue = $optDataForPlotHsh_ref->{$dataType}{$opt_minTransfragLen}{$opt_minTransfragPosCov}{$opt_minTransfragCovPerNt}{$opt_maxTransfragGapTorlerance}{$opt_splitWinSize}{$opt_splitFold};
							$extraStatment = "optimalPoint <- data.frame(maxTransfragGapTorlerance=$opt_maxTransfragGapTorlerance, $dataType=$optimalValue, splitFold=$opt_splitFold);";
							$extraArg = " + geom_point(data = optimalPoint, colour=\"green\", size = 20, alpha = 0.3)";
						}
						
						if ($dataType ne 'score') {
							$extraArg .= " + geom_hline(aes(yintercept=$dataTypeStatHsh_ref->{$dataType}{'mean'}), colour=\"black\", size = 1, linetype = \"dashed\")";
						}
					
						my $nameTag = $transfragLenCovTag;
						my $dataPath = $ggplotDirHsh_ref->{'dat'}."/$dirTag/$nameTag.dat";
						my $pdfPath = $ggplotDirHsh_ref->{'pdf'}."/$dirTag/$nameTag.pdf";
						my $RScriptPath = $ggplotDirHsh_ref->{'R'}."/$dirTag/$nameTag.R";
						my $logPath = $ggplotDirHsh_ref->{'log'}."/$dirTag/$nameTag.log";
						my $factorA = 'maxTransfragGapTorlerance';
						my $factorB = 'splitWinSize';
						my $factorC = 'splitFold';
						my $YAxis = $dataType;
						my $height = 5;
						my $width = 15;
						&ggplotXYLineThreeFactors($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $factorA, $factorB, $factorC, $YAxis, $extraStatment, $extraArg, $height, $width);#->1039
					}
				}
			}
		}
	}
	
	return ();
}
sub printBothFHAndStdout {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: currentTime|818
#	appearInSub: calculateOptimalTransfragScore|276
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_optimizeParameters|185
#	input: $FH, $message, $numTrailingSpace
#	output: 
#	toCall: &printBothFHAndStdout($message, $numTrailingSpace, $FH);
#	calledInLine: 361, 362, 363, 365
#....................................................................................................................................................#
	
	my ($message, $numTrailingSpace, $FH) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);

	print "[".&currentTime()."] ".$message.$trailingSpaces."\n";#->818
	print {$FH} "[".&currentTime()."] ".$message."\n";#->818

	return ();
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|818
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|65, 9_finishingTasks|218
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 71, 224
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;

	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->818
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->818
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->818
		print "=========================================================================\n\n";
	}
}
sub printGFF_oneRNAPerGene_chooseStrnd_filter {
#....................................................................................................................................................#
#	subroutineCategory: gff
#	dependOnSub: >none
#	appearInSub: outputFinalTransfragGff|1311
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 8_outputFiles|207
#	input: $filter, $geneInfoHsh_ref, $gffGeneLineOnly, $outGFFPath, $strndInWord
#	output: none
#	toCall: &printGFF_oneRNAPerGene_chooseStrnd_filter($geneInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, $filter);
#	calledInLine: 1336
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, $filter) = @_;

	#------gffGeneLineOnly: print only the gene line, which is perfect for trnsfrgs
	#------strndInWord: plus, minus or both
		
	$outGFFPath =~ s/\/+/\//g;
	
	my $strndOut = '+';
	$strndOut = '-' if $strndInWord eq 'minus';
	
	open (GFFOUT, ">", "$outGFFPath");
	print GFFOUT "##gff-version\t3\n";

	foreach my $geneID (sort {$a cmp $b} keys %{$geneInfoHsh_ref}) {
		
		if ($filter) {
			if ($geneInfoHsh_ref->{$geneID}{$filter}) {
				next if $geneInfoHsh_ref->{$geneID}{$filter} eq 'no';
			}
		}
		
		my $strnd = $geneInfoHsh_ref->{$geneID}{"strnd"};

		next if $strndOut ne $strnd and $strndInWord ne 'both';
		
		my $cntg = $geneInfoHsh_ref->{$geneID}{"cntg"};
		my $ctgry = $geneInfoHsh_ref->{$geneID}{"ctgry"};
		my $description = $geneInfoHsh_ref->{$geneID}{"description"};
		
		my ($geneStart, $geneEnd) = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{'geneRng'}};
		print GFFOUT join "", (join "\t", ($cntg, 'BCP', 'gene', $geneStart, $geneEnd, ".", $strnd, ".", "ID=$geneID;Name=$geneID;description=$description")), "\n";

		if ($gffGeneLineOnly eq 'no') {#-----will print the RNA and exon lines also, aim to avoid display annoyance on IGV
			my $RNAID = $geneInfoHsh_ref->{$geneID}{"RNAID"};
			my ($RNAStart, $RNAEnd) = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{'RNARng'}};
			print GFFOUT join "", (join "\t", ($cntg, 'BCP', $ctgry, $RNAStart, $RNAEnd, ".", $strnd, ".", "ID=$RNAID;Name=$geneID;Parent=$geneID;description=$description;")), "\n";
			
			foreach my $rngType ('exonRng', 'CDSRng') {
				my $gffRng = 'exon';
				$gffRng = 'CDS' if $rngType eq 'CDSRng';
				
				if ($geneInfoHsh_ref->{$geneID}{$rngType}) {
					my @rngAry = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{$rngType}};
					my $num = 1;
					for (my $i=0; $i < $#rngAry; $i += 2) {
						$num++;
						my ($start, $end) = ($rngAry[$i], $rngAry[$i+1]);
						my $ID = "$rngType\_$num\_$RNAID";
						print GFFOUT join "", (join "\t", ($cntg, "BCP", $gffRng, $start, $end, ".", $strnd, ".", "ID=$ID;Name=$ID;Parent=$RNAID;description=.;")), "\n";
					}
				}
			}
		}
	}
	close GFFOUT;
}
sub readGFF_oneRNAPerGene {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: currentTime|818
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|163
#	secondaryAppearInSection: >none
#	input: $gffPath
#	output: $geneInfoHsh_ref
#	toCall: my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);
#	calledInLine: 172
#....................................................................................................................................................#

	my ($gffPath) = @_;

	my $geneInfoHsh_ref = {};

	#---read the gff
	my $geneByRNAHsh_ref = {};

	open (GFF, $gffPath);
	print "[".&currentTime()."] Reading: $gffPath\n";#->818
	while (my $theLine = <GFF>) {

		chomp $theLine;

		last if $theLine =~ m/^##FASTA/;

		if ($theLine !~ m/^\#|^\@/ and $theLine !~ m/\tsupercontig\t/) {

			my ($seq, undef, $geneCategory, $featureStart, $featureEnd, undef, $geneStrd, undef, $dscrptns) = split (/\t/, $theLine);

			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));

			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene

				my $geneID = $unqID;

				$geneInfoHsh_ref->{$geneID}{'strnd'} = $geneStrd;
				$geneInfoHsh_ref->{$geneID}{'cntg'} = $seq;
				$geneInfoHsh_ref->{$geneID}{'description'} = uri_unescape($geneName);
				$geneInfoHsh_ref->{$geneID}{'description'} =~ s/\+/ /g;
				@{$geneInfoHsh_ref->{$geneID}{'geneRng'}} = ($featureStart, $featureEnd);

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes

				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'CDSRng'}}, ($featureStart, $featureEnd);

			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'exonRng'}}, ($featureStart, $featureEnd);

			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh_ref->{$RNAID} = $geneID;
				$geneInfoHsh_ref->{$geneID}{'ctgry'} = $geneCategory;
				@{$geneInfoHsh_ref->{$geneID}{'RNARng'}} = ($featureStart, $featureEnd);
				$geneInfoHsh_ref->{$geneID}{'RNAID'} = $RNAID;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close GFF;

	#---get the UTR if any
	my $minUTRLength = 10;
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {

		if (exists $geneInfoHsh_ref->{$geneID}{'CDSRng'}) {
			my $exonMin = min(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $exonMax = max(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $CDSMin = min(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});
			my $CDSMax = max(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});

			if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			} else {
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			}
		}
	}

	return ($geneInfoHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|65
#	secondaryAppearInSection: >none
#	input: none
#	output: $gffPath, $outDir, $queryCovPlsIndexPath
#	toCall: my ($queryCovPlsIndexPath, $gffPath, $outDir) = &readParameters();
#	calledInLine: 72
#....................................................................................................................................................#

	my ($queryCovPlsIndexPath, $gffPath, $outDir);
	my $dirPath = dirname(rel2abs($0));

	$outDir = "$dirPath/transfragDiscoverer/";

	GetOptions 	("queryCovPlsIndexPath=s" => \$queryCovPlsIndexPath,
				 "gffPath=s"  => \$gffPath,
				 "outDir:s"  => \$outDir)

	or die		("Error in command line arguments\n");

	#---check file
	foreach my $fileToCheck ($queryCovPlsIndexPath, $gffPath) {
		die "Can't read $fileToCheck" if (not -s $fileToCheck and not -s "$fileToCheck.gz");
	}

	system "mkdir -p -m 777 $outDir/";

	return($queryCovPlsIndexPath, $gffPath, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|818
#	appearInSub: checkGeneInfo|372, checkOverlapAndProximity|398, checkmRNAProximity|633, checkmRNATrnsfrgOverlap|664, classifymRNABasedOnProximity|770, getmRNASubset|1008, identifyTransfrag|1078, optimizeTransfrgIdentificationParameters|1254, plotTransfrgIdentificationParameters|1424
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_processInputData|163, 6_optimizeParameters|185
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 386, 394, 443, 457, 469, 649, 652, 678, 684, 787, 811, 1034, 1105, 1217, 1249, 1286, 1289, 1449
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);

	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->818

	return ();
}
sub scanDramaticCoverageChange {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: identifyTransfrag|1078
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $minCovPerNtInWindow, $splitFold, $splitStepSize, $splitWinSize, $strnd, $tmpTrnsfrgCovHsh_ref, $tmpTrnsfrgPosHsh_ref
#	output: $splitTrnsfrgRngHsh_ref
#	toCall: my ($splitTrnsfrgRngHsh_ref) = &scanDramaticCoverageChange($tmpTrnsfrgPosHsh_ref, $tmpTrnsfrgCovHsh_ref, $strnd, $splitWinSize, $splitFold, $minCovPerNtInWindow, $splitStepSize);
#	calledInLine: 1184
#....................................................................................................................................................#

	my ($tmpTrnsfrgPosHsh_ref, $tmpTrnsfrgCovHsh_ref, $strnd, $splitWinSize, $splitFold, $minCovPerNtInWindow, $splitStepSize) = @_;

	my $splitTrnsfrgRngHsh_ref = {};

	if ($splitWinSize <= 0) {#----split trnsfrg turned off

		my $start = 0;
		my $end = $#{$tmpTrnsfrgCovHsh_ref};
		$splitTrnsfrgRngHsh_ref->{$start} = $end;

	} else {#-----split the trnsfrg

		my $halfWinSize = int ($splitWinSize/2);
		my $trnsfrgLen = @{$tmpTrnsfrgPosHsh_ref->{$strnd}};
		my $initialBreakPtAry_ref = ();
		my $finalBreakPtAry_ref = ();
		my $splitNum = 0;

		#---sliding window in $splitStepSize
		for (my $i = $halfWinSize; $i <= $trnsfrgLen-$halfWinSize-1; $i += $splitStepSize) {

			my $tmpFoldChange = 0;

			my $leftWinSum = sum(@{$tmpTrnsfrgCovHsh_ref->{$strnd}}[$i-$halfWinSize..$i]);
			my $rightWinSum = sum(@{$tmpTrnsfrgCovHsh_ref->{$strnd}}[$i..$i+$halfWinSize]);

			#---do comparison only if the covPerNt if great than minCovPerNtInWindow
			if ((($leftWinSum+$rightWinSum)/$splitWinSize) >= $minCovPerNtInWindow) {
			
				if ($leftWinSum > $rightWinSum) {
					eval {$tmpFoldChange = $leftWinSum/$rightWinSum;};
				} elsif ($rightWinSum > $leftWinSum) {
					eval {$tmpFoldChange = $rightWinSum/$leftWinSum;};
				}

				if ($tmpFoldChange >= $splitFold) {
					push @{$initialBreakPtAry_ref}, $i;
					$splitNum++;
				}
			}
		}

		#---group the continous breakpoints into clusters
		push @{$finalBreakPtAry_ref}, 0;#----push the first position

		if ($splitNum > 0) {#---there are breakpoints

			my $collapsedBreakPtHsh_ref = {};
			my $clusterNum = 0;
			push @{$collapsedBreakPtHsh_ref->{$clusterNum}}, $initialBreakPtAry_ref->[0];

			#---more than 1 breakpoint
			if (@{$initialBreakPtAry_ref} > 1) {
				foreach my $i (1..$#{$initialBreakPtAry_ref}) {
					$clusterNum++ if (($initialBreakPtAry_ref->[$i] - $initialBreakPtAry_ref->[$i-1]) > $splitWinSize); #---if distance between breakpoints > ($splitWinSize) treated as another cluster;
					push @{$collapsedBreakPtHsh_ref->{$clusterNum}}, $initialBreakPtAry_ref->[$i];
				}
			}

			#---collapse the breakpoint clusters
			foreach my $clusterNum (sort {$a<=>$b} keys %{$collapsedBreakPtHsh_ref}) {
				my $clusterMidPt = int(sum(@{$collapsedBreakPtHsh_ref->{$clusterNum}})/@{$collapsedBreakPtHsh_ref->{$clusterNum}});
				push @{$finalBreakPtAry_ref}, $clusterMidPt;
				push @{$finalBreakPtAry_ref}, $clusterMidPt+1;
			}
		}

		push @{$finalBreakPtAry_ref}, $#{$tmpTrnsfrgPosHsh_ref->{$strnd}};#----push the last position

		@{$finalBreakPtAry_ref} = sort {$a <=> $b} @{$finalBreakPtAry_ref}; #--- the @finalRngAry looks like (0, $clusterMidPt1, $clusterMidPt1+1, $clusterMidPt2, $clusterMidPt2+1, $transfragEnd) if there are two collapsed break points;

		#---convert into hash to report
		foreach (my $i=0; $i < $#{$finalBreakPtAry_ref}; $i += 2) {#---(0, 2, 4, 8.....)
			my $start = ${$finalBreakPtAry_ref}[$i];
			my $end = ${$finalBreakPtAry_ref}[$i+1];
			$splitTrnsfrgRngHsh_ref->{$start} = $end;
		}
	}

	return ($splitTrnsfrgRngHsh_ref);

}

exit;
