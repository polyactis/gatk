package org.broadinstitute.sting.gatk.walkers.yh;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;
import java.util.regex.Pattern;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.yh.MergeVCFReplicateHaplotypes.Diplotype;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

public class CalculateConcordanceAmongReplicates extends MergeVCFReplicateHaplotypes {
	private @Argument(fullName = "concordanceStatFname", shortName = "concordanceStatFname", 
			doc = "output file to contain concordance statistics",
			required = true)
	String concordanceStatFname;
	//private @Output(doc = "output file to contain concordance results", required = true)
	private FileWriter concordanceStatWriter;
	
	
	Map<String, List<String>> trueSampleID2SampleIDList = new HashMap<String, List<String>>();

	public void initialize() {
		/*
		 * 2011-11-11
		 * 	things to initialize before the mapper is called
		 */
		
		Pattern myPattern = Pattern.compile(replicateIndividualTag);
		
		List<String> rodNames = Arrays.asList(variantCollection.variants.getName());
		Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
		
		//Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), null);

		TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
		Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);
		
		Iterator<String> itr = vcfSamples.iterator();
		while (itr.hasNext()){
			String sampleID = (String) itr.next();
			String trueSampleID = sampleID.split(myPattern.toString())[0];
			trueSampleIDSet.add(trueSampleID);
			sampleID2trueSampleID.put(sampleID, trueSampleID);
			if (!trueSampleID2Count.containsKey(trueSampleID)){
				trueSampleID2Count.put(trueSampleID, 0);
				trueSampleID2SampleIDList.put(trueSampleID, new ArrayList<String>());
			}
			trueSampleID2Count.put(trueSampleID, trueSampleID2Count.get(trueSampleID)+1);
			trueSampleID2SampleIDList.get(trueSampleID).add(sampleID);
		}
		
		
		//2012.7.23
		if (concordanceStatFname != null) {
			try {
				concordanceStatWriter = new FileWriter(concordanceStatFname);
			} catch (IOException e) {
				e.printStackTrace();
				System.err.println("Error: " + e.getMessage());
			}
			//2012.7.23 output the majority vote support & replicate count/family size for each replicate
			Vector<String> header = new Vector<String>();
			header.add("replicate1ID");
			header.add("replicate2ID");
			header.add("noOfMatches");
			header.add("noOfLoci");
			header.add("concordance");
			header.add("noOfMatchesHomo");
			header.add("noOfLociHomo");
			header.add("concordanceHomo");
			outputStringVector(concordanceStatWriter, header);
		
		}
	}
	
	
	public void onTraversalDone(Integer result) {
		logger.info(result + " loci to be processed.");
		logger.info(sampleID2Diplotype.size() + " samples.");
		logger.info("");
		
		int i, j, k;
		String replicate1ID, replicate2ID;
		byte[] replicate1Genotype, replicate2Genotype;
		Diplotype replicate1Diplotype, replicate2Diplotype;
		Integer _noOfMatches, _noOfLoci, _noOfMatchesHomo, _noOfLociHomo;
		Double concordance, concordanceHomo;
		int noOfReplicatePairs = 0;
		int noOfUniqueSamplesWithReplicates = 0;
		for (Map.Entry<String, List<String>> item: trueSampleID2SampleIDList.entrySet()){
			// loop through each sample (column or individual)
			String trueSampleID = item.getKey();
			List<String> sampleIDList = item.getValue();
			if (sampleIDList.size()>=2){
				noOfUniqueSamplesWithReplicates++;
			}
			for (i=0;i<sampleIDList.size()-1;i++){
				replicate1ID = sampleIDList.get(i);
				replicate1Diplotype = sampleID2Diplotype.get(replicate1ID);
				for (j=i+1; j<sampleIDList.size();j++){
					noOfReplicatePairs ++;
					replicate2ID = sampleIDList.get(j);
					replicate2Diplotype = sampleID2Diplotype.get(replicate2ID);
					_noOfMatches =0;
					_noOfLoci = 0;
					_noOfMatchesHomo = 0;	//homozygotes at both replicates
					_noOfLociHomo = 0;	//homozygotes at both replicates
					for (k=0; k<replicate2Diplotype.noOfLoci; k++){
						replicate1Genotype = replicate1Diplotype.getGenotypeAtOnePosition(k+1, true);
						replicate2Genotype = replicate2Diplotype.getGenotypeAtOnePosition(k+1, true);
						if (replicate1Genotype[0]==replicate2Genotype[0] && replicate1Genotype[1]==replicate2Genotype[1]){
							_noOfMatches++;
							if (replicate1Genotype[0]==replicate1Genotype[1] && replicate2Genotype[0]==replicate2Genotype[1]){
								_noOfMatchesHomo++;
							}
						}
						_noOfLoci++;
						if (replicate1Genotype[0]==replicate1Genotype[1] && replicate2Genotype[0]==replicate2Genotype[1]){
							_noOfLociHomo++;
						}
						
					}
					concordance =  _noOfMatches.doubleValue()/_noOfLoci.doubleValue();
					concordanceHomo = _noOfMatchesHomo.doubleValue()/_noOfLociHomo.doubleValue();
					
					//output the stat
					Vector<String> statDataRow  = new Vector<String>();
					statDataRow.add(replicate1ID);
					statDataRow.add(replicate2ID);
					statDataRow.add(_noOfMatches.toString());
					statDataRow.add(_noOfLoci.toString());	
					statDataRow.add(concordance.toString());
					statDataRow.add(_noOfMatchesHomo.toString());
					statDataRow.add(_noOfLociHomo.toString());
					statDataRow.add(concordanceHomo.toString());
					
					outputStringVector(concordanceStatWriter, statDataRow);	
				}
			}
		}
		
		//2012.6.12 release memory
		sampleID2Diplotype.clear();
		
		//2012.7.24 close the file, otherwise it might end up with some missing end (buffer not flushed to disk)
		closeFileWriter(concordanceStatWriter);
		
		// 2012.7.24 don't try to close vcfWriter (it'll close itself somewhere after)
		logger.info(result + " loci processed.");
		logger.info(noOfUniqueSamplesWithReplicates + " unique samples with replicates.");
		logger.info(noOfReplicatePairs + " replicate pairs.");
		
		logger.info("Number of total genotypes: " + no_of_total_genotypes);
		logger.info("Number of outputted genotypes: "
				+ no_of_included_genotypes);
		logger.info("Number of non-bi-allelic loci after GQ and Depth filter: "
				+ no_of_nonBiAllelicLoci);
	}
}
