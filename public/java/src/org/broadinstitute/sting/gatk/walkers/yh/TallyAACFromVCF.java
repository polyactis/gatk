/*
 * 2011.11.17
 * 	Yu Huang
 */
package org.broadinstitute.sting.gatk.walkers.yh;

import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

/**
 **	2011-11-15 a program that count the number of loci with particular alternative allele count.
 *		If AC/AC1 is not available in INFO field, it'll calculate it from the genotypes.
 *	Note:
 *	
 *	Examples:
 *	java  -Xmx2g -jar /home/crocea/script/gatk/dist/GenomeAnalysisTK.jar 
 *		-R  /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
 *		-T TallyAACFromVCF --variant samtools/Contig0.vcf.gz
 *		-o tmp/SAMtoolsConitg0_TallyAAC.tsv
 *
 */
public class TallyAACFromVCF  extends RodWalker<Integer, Integer>{
	@ArgumentCollection protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();
	
	@Output(doc="output filename to which tally stats should be written", required=true)
	public String outputFname;
	
	Map<Integer, Integer> AAC2Count = new HashMap<Integer, Integer>();
	Integer no_of_total_genotypes = 0;
	Integer no_of_included_genotypes = 0;
	Integer no_of_retained_loci = 0;
	
	String[] outputHeader = new String[2];
	FileWriter writer;
	FileOutputStream outStream;
	public void initialize() {
		/*
		 * 2011-11-11
		 * 	things to initialize before the mapper is called
		 */
		try {
			writer = new FileWriter(outputFname);
			outputHeader[0] = "AAC";
			outputHeader[1] = "NumberOfLoci";
			for(int i = 0; i < outputHeader.length; i++)
			{
				writer.append(outputHeader[i]);
				if(i < outputHeader.length - 1)
					writer.append("\t");
				else
					writer.append('\n');
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.err.println("Error: " + e.getMessage());
		}
	}
	
	@Override
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		// TODO Auto-generated method stub
		if (tracker==null){
			return 0;
		}
		Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());
		/*
		 * 2011-11-11 line from the example in the 20linesaver slides doesn't work here
		 */
		//Collection<VariantContext> vcs = tracker.getVariantContexts(ref, "variant", null, context.getLocation(), true, true);

		if (vcs == null || vcs.size() == 0) {
			return 0;
		}
		
		
		for (VariantContext vc: vcs){
			int AC = -1;
			if (vc.hasAttribute("AC")
					&& vc.getAttribute("AC") instanceof Integer) {
				AC = vc.getAttributeAsInt("AC", 0);
			}
			else if (vc.hasAttribute("AC1")
					&& vc.getAttribute("AC1") instanceof Integer) {
				AC = vc.getAttributeAsInt("AC1", 0);
			}
			else if (vc.isVariant()) {
				for (Allele allele : vc.getAlternateAlleles())
					AC = Math.max(AC, vc.getCalledChrCount(allele));
			} else
				// by default, the site is considered monomorphic
				AC = 0;
			if (AAC2Count.containsKey(AC)){
				AAC2Count.put(AC, AAC2Count.get(AC)+1);
			}
			else{
				AAC2Count.put(AC, 1);
			}
			no_of_total_genotypes += vc.getNSamples();
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
		return value+sum;
	}
	
	public void onTraversalDone(Integer result) {
		try {
			for (Map.Entry<Integer, Integer> item : AAC2Count.entrySet()) {
					//System.err.println("AAC: " + item.getKey());
					//System.err.println("NumberOfLoci: " + item.getValue());
					writer.write(item.getKey().toString());
					writer.write("\t");
					writer.write(item.getValue().toString());
					writer.write("\n");
				
				
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.err.println("Error: " + e.getMessage());
		}
		logger.info(result + " loci processed.");
		logger.info("Number of total genotypes: " + no_of_total_genotypes);
	}
	
}
