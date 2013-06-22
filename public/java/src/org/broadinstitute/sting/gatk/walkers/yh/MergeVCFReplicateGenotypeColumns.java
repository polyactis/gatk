/*
 * 2012.3.30
 * 	a RODwalker to merge unphased calls that are from the same individual.
 * 		part of the AlignmentToTrioCallPipeline.py workflow.
 * 
 *  * Note: 1. Samples from VCF that are not in depthFile would be ignored in out.
 * 2. AC/AF/DP on the INFO column will be updated but QUAL will not.
 * 
 * Examples (--onlyKeepBiAllelicSNP) is optional: java -Xmx2g -jar
 * /home/crocea/script/gatk/dist/GenomeAnalysisTK.jar -T MergeVCFReplicateGenotypeColumns
 *  -R /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
 *  --variant samtools/Contig0.vcf.gz -o /tmp/contig0_afterMerge.vcf --onlyKeepBiAllelicSNP
 *  --replicateIndividualTag copy
 */
package org.broadinstitute.sting.gatk.walkers.yh;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

public class MergeVCFReplicateGenotypeColumns extends MergeVCFReplicateHaplotypes{
	//@ArgumentCollection
	//protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();
	
	// key is sample ID, value is a list of minDepth, maxDepth, [minGQ|-1.0].
	private TreeSet<String> samples = new TreeSet<String>();

	
	// 2012.7.30 no initialize() , inherited from MergeVCFReplicateHaplotypes
	
	public class ReplicateGenotypeHolder{
		private Map<String, ArrayList<Genotype>> GT2GenotypeList;
		private String sampleID;
		Integer maxSupport=0;
		String GT_withMaxSupport = null;
		Genotype genotypeWithMaxSupport = null;
		public ReplicateGenotypeHolder(String sampleID){
			this.sampleID = sampleID;
			this.GT2GenotypeList = new HashMap<String, ArrayList<Genotype>>();
		}
		
		public void addGenotype(Genotype genotype){
			String gt = genotype.getGenotypeString();
			//String alleleA = genotype.getAllele(0).getBaseString();
			
			//Allele originalAlleleB = (g.getAlleles().size() == 2) ? g.getAllele(1) : g.getAllele(0); // hack to deal with no-call genotypes
			//Genotype imputedGenotype = new Genotype(g.getSampleName(), alleles, genotypeQuality, filters,originalAttributes , genotypeIsPhased);
			if (!this.GT2GenotypeList.containsKey(gt)){
				this.GT2GenotypeList.put(gt, new ArrayList<Genotype>());
			}
			this.GT2GenotypeList.get(gt).add(genotype);
		}
		public Genotype getGenotypeWithMaxSupport(){
			for (Map.Entry<String, ArrayList<Genotype>> item: GT2GenotypeList.entrySet() ){
				Integer support = item.getValue().size();
				if (support>maxSupport){
					maxSupport = support;
					GT_withMaxSupport = item.getKey();
				}
			}
			if (maxSupport>0){
				genotypeWithMaxSupport =  GT2GenotypeList.get(GT_withMaxSupport).get(0);	//return the 1st genotype
				
			}
			return genotypeWithMaxSupport;
		}
	}
	
	@Override
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		if (tracker==null){
			return 0;
		}
		
		Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());
		
		if (vcs == null || vcs.size() == 0) {
			return 0;
		}
		
		for (VariantContext vc: vcs){
			Map <String, ReplicateGenotypeHolder> trueSampleID2ReplicateGenotypeHolder = new HashMap<String, ReplicateGenotypeHolder>();
			GenotypesContext thisGC = vc.getGenotypes();
			for (Genotype genotype: thisGC) {
				no_of_total_genotypes ++;
				String sampleID = genotype.getSampleName();
				String trueSampleID = sampleID2trueSampleID.get(sampleID);
				if (!trueSampleID2ReplicateGenotypeHolder.containsKey(trueSampleID)){
					trueSampleID2ReplicateGenotypeHolder.put(trueSampleID, new ReplicateGenotypeHolder(trueSampleID));
				}
				trueSampleID2ReplicateGenotypeHolder.get(trueSampleID).addGenotype(genotype);
				
			}
			
			Set<String> selectedSamples = new HashSet<String>(0);
			GenotypesContext newGC = GenotypesContext.create();
			//newGC.clear();
			for (Map.Entry<String, ReplicateGenotypeHolder> item: trueSampleID2ReplicateGenotypeHolder.entrySet()){
				Genotype genotypeWithMaxSupport = item.getValue().getGenotypeWithMaxSupport();
				selectedSamples.add(genotypeWithMaxSupport.getSampleName());
				// change its sampleID
				genotypeWithMaxSupport = Genotype.modifyName(genotypeWithMaxSupport, item.getKey());
				newGC.add(genotypeWithMaxSupport);
			}
			
			VariantContext sub = vc.subContextFromSamples(selectedSamples, vc.getAlleles());
			VariantContextBuilder builder = new VariantContextBuilder(sub);
			builder.genotypes(newGC);
			int depth = 0;
			for (Genotype g : newGC) {
				if (g.isNotFiltered() && g.isCalled()) {
					Integer dp = getGenotypeDepth(g);
					depth += dp;
				}
			}

			if (KEEP_ORIGINAL_CHR_COUNTS) {
				if (sub.hasAttribute(VCFConstants.ALLELE_COUNT_KEY))
					builder.attribute("AC_Orig",
							sub.getAttribute(VCFConstants.ALLELE_COUNT_KEY));
				if (sub.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY))
					builder.attribute("AF_Orig",
							sub.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
				if (sub.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY))
					builder.attribute("AN_Orig",
							sub.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
			}

			VariantContextUtils.calculateChromosomeCounts(builder, false);
			builder.attribute("DP", depth);

			VariantContext newVC = builder.make();
			
			
			
			if (onlyKeepBiAllelicSNP==true){
				//2011.11.28
				if (newVC.getAlleles().size()!=2){
					no_of_nonBiAllelicLoci ++;
					continue;
				}
			}
			if (newVC.getGenotypes().size()>0){	//output only when not all are missing
				no_of_retained_loci ++;
				if (vcfWriter!=null){
					vcfWriter.add(newVC);
				}
			}
		}
		return 1;
	}

	@Override
	public Integer reduceInit() {
		return 0;
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
		return value+sum;
	}
	/**
	 * 2012.3.30 copied from FilterVCFByDepth.java
	 * 
	 * @param genotype
	 * @return
	 */
	public Integer getGenotypeDepth(Genotype genotype) {
		Integer genotypeDepth;
		String dp = (String) genotype.getAttribute("DP");
		if (dp != null && !dp.equals(VCFConstants.MISSING_DEPTH_v3)
				&& !dp.equals(VCFConstants.MISSING_VALUE_v4)) {
			genotypeDepth = Integer.valueOf(dp);
		} else {
			genotypeDepth = 0;
		}
		return genotypeDepth;
	}
	
	public void onTraversalDone(Integer result) {
		logger.info(result + " loci processed.");
		logger.info(no_of_retained_loci + " loci retained.");
		logger.info("Number of total genotypes: " + no_of_total_genotypes);
		logger.info("Number of outputted genotypes: "
				+ no_of_included_genotypes);
		logger.info("Number of non-bi-allelic loci after GQ and Depth filter: "
				+ no_of_nonBiAllelicLoci);
	}

}
