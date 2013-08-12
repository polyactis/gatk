/*
 * 2013.08.12 a program that filters loci in VCF file based on the heterozygous fraction of each locus
 *		java -jar dist/GenomeAnalysisTK.jar -T FilterVCFBySiteHetFraction
 *		--variant /Network/Data/vervet/db/genotype_file/method_101/89279_VCF_88057_VCF_Scaffold3.filterByMaxSNPMissingRate.recode.vcf.gz
 *		--maxHetFraction 0.6 --reference_sequence /Network/Data/vervet/db/individual_sequence/3280_vervet_ref_6.0.3.fasta
 *		--out /tmp/Scaffold3.maxHetFraction0.6.vcf
 *
 */
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

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

public class FilterVCFBySiteHetFraction  extends RodWalker<Integer, Integer> {
	@ArgumentCollection
	protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

	@Output(doc = "File to which variants should be written", required = false)
	public VCFWriter vcfWriter;

	@Argument(fullName = "maxHetFraction", shortName = "mHF", doc = "maximum fraction of heterozygous calls in one site.",
			required = false)
	double maxHetFraction = 0.6;
	
	Long no_of_total_loci = (long) 0;
	Long no_of_retained_loci = (long) 0;
	public void initialize() {
		
		// 2012.5.1 get all samples from the vcf input file and headers
		List<String> rodNames = Arrays.asList(variantCollection.variants.getName());
		Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
		
		Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);
		TreeSet<String> vcfSamples = new TreeSet<String>(
				SampleUtils.getSampleList(vcfRods,
						VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
		
		if (vcfWriter != null) {
			// 2011-11-11 borrow the header from the input VCF
			vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));
			// VCFUtils.getHeaderFields(this.getToolkit(),
			// SampleUtils.getUniqueSamplesFromRods(this.getToolkit()) returns
			// empty.
		}
	}
	
	@Override
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		// TODO Auto-generated method stub
		if (tracker == null) {
			return 0;
		}
		Collection<VariantContext> vcs = tracker.getValues(
				variantCollection.variants, context.getLocation());
		
		if (vcs == null || vcs.size() == 0) {
			return 0;
		}
		float no_of_non_missing_genotypes = 0;
		int no_of_hets = 0;
		for (VariantContext vc : vcs) {
			no_of_total_loci += 1;	//vc.getCalledChrCount();
				
				
			GenotypesContext thisGC = vc.getGenotypes();
			// Map <String, Genotype> sampleID2genotype = vc.getGenotypes();
			ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
			no_of_non_missing_genotypes = 0;
			no_of_hets = 0;
			for (Genotype genotype : thisGC) {
				if (genotype.isAvailable() && genotype.isCalled()){
					no_of_non_missing_genotypes += 1;
					if (genotype.isHet()){
						no_of_hets++;
					}
				}
			}
			float hetFraction = no_of_hets/no_of_non_missing_genotypes;
			
				
			if (hetFraction<=maxHetFraction) { 
				no_of_retained_loci++;
				//sub = subsetRecord(vc, selectedSamples);
				if (vcfWriter != null) {
					vcfWriter.add(vc);
				}
			}
			
		}
		return 1;
	}
	
	@Override
	public Integer reduceInit() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
		return value + sum;
	}

	public void onTraversalDone(Integer result) {
		logger.info(result + " loci processed.");
		logger.info("Number of total loci: " + no_of_total_loci);
		logger.info("Number of loci retained: " + no_of_retained_loci );
	}
}
