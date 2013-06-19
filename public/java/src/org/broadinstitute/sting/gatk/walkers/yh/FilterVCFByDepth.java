/*
 * 2011-11-11
 * 	a RODwalker to mask calls as missing if it's below or above a range.
 */
package org.broadinstitute.sting.gatk.walkers.yh;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
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
import java.util.Vector;

import org.broad.tribble.readers.LineReader;
import org.broad.tribble.util.ParsingUtils;
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
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

/**
 ** 2011-11-15 a program that mark calls as missing based on
 * minDepth/maxDepth/minGQ from depthFile. this program stole quite a bit from
 * SelectVariants Format of tab-delimited depthFile: sampleID/readGroup minDepth
 * maxDepth minGQ
 * 
 * Note: 1. Samples from VCF that are not in depthFile would be ignored in out.
 * 2. AC/AF/DP on the INFO column will be updated but QUAL will not.
 * 
 * Examples (--onlyKeepBiAllelicSNP) is optional: java -Xmx2g -jar
 * /home/crocea/script/gatk/dist/GenomeAnalysisTK.jar -T FilterVCFByDepth
 *  -R /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
 *  --variant samtools/Contig0.vcf.gz -depthFname /tmp/coverage.tsv
 *  -o tmp/SAMtoolsConitg0_afterDepthFilter.vcf --onlyKeepBiAllelicSNP
 *  -ssFname /tmp/siteStatContig0.tsv
 * 
 */
public class FilterVCFByDepth extends RodWalker<Integer, Integer> {
	@ArgumentCollection
	protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

	@Output(doc = "File to which variants should be written", required = false)
	public VCFWriter vcfWriter;

	@Argument(fullName = "sampleMinMaxDepthFname", shortName = "depthFname", doc = "tab-delimited file, sample-id, minDepth, maxDepth", required = true)
	String sampleMinMaxDepthFname;

	@Argument(fullName = "siteStatFname", shortName = "ssFname", doc = "this file (if given) would contain chr, pos, missing-rate, AC, AF. tab-delimited.", required = false)
	String siteStatFname;
	FileWriter siteStatWriter;
	Vector<String> siteStatOutputHeader = new Vector<String>();

	// 2011.11.28
	@Argument(fullName = "onlyKeepBiAllelicSNP", doc = "set this to true to remove SNPs with >=3 alleles", required = false)
	Boolean onlyKeepBiAllelicSNP = false;

	Map<String, List<Float>> sampleID2minMaxDepth = new HashMap<String, List<Float>>();
	Float totalMinDepth = (float) 0;	//2012.5.1 sum of minDepth from the sampleID2minMaxDepth and the sample has to be in VCF headers.
	Float totalMaxDepth = (float) 0;	//2012.5.1
	
	// key is sample ID, value is a list of minDepth, maxDepth, [minGQ|-1.0].
	Long no_of_total_genotypes = (long) 0;
	Long no_of_included_genotypes = (long) 0;
	Long no_of_retained_loci = (long) 0;
	Long no_of_nonBiAllelicLoci = (long) 0;
	private TreeSet<String> samples = new TreeSet<String>();

	// 2012.3.30 do keep the original AC, AF, AN
	private boolean KEEP_ORIGINAL_CHR_COUNTS = true;

	public void initialize() {
		/*
		 * 2011-11-11 things to initialize before the mapper is called
		 */
		// 2011.11.28 report this variable
		logger.info("onlyKeepBiAllelicSNP= " + onlyKeepBiAllelicSNP + "\n");
		
		// 2012.5.1 get all samples from the vcf input file and headers
		List<String> rodNames = Arrays.asList(variantCollection.variants
				.getName());
		Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(
				getToolkit(), rodNames);

		// Map<String, VCFHeader> vcfRods =
		// VCFUtils.getVCFHeadersFromRods(getToolkit(), null);

		TreeSet<String> vcfSamples = new TreeSet<String>(
				SampleUtils.getSampleList(vcfRods,
						VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
		Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(
				vcfRods.values(), logger);

		Iterator<String> itr = vcfSamples.iterator();
		while (itr.hasNext()) {
			String sampleID = (String) itr.next();
			//if (sampleID2minMaxDepth.containsKey(sampleID)) {
			//}
			samples.add(sampleID);
		}
		logger.info(samples.size() + " samples in the input VCF.");
		
		
		File sampleMinMaxDepthFile = new File(sampleMinMaxDepthFname);
		FileInputStream sampleMinMaxDepthInStream;
		int noOfSamplesWithDepth = 0;
		try {
			sampleMinMaxDepthInStream = new FileInputStream(
					sampleMinMaxDepthFname);

			if (!sampleMinMaxDepthFile.canRead()) {
				logger.error(sampleMinMaxDepthFile.getAbsolutePath()
						+ " is not readable.");
				System.err.println(sampleMinMaxDepthFile.getAbsolutePath()
						+ " is not readable.\n");
				System.exit(0);
			}
			String line;
			Integer lineNo = 0;
			BufferedReader br = new BufferedReader(new InputStreamReader(
					sampleMinMaxDepthInStream));
			line = br.readLine(); // skip the header
			while ((line = br.readLine()) != null) {
				String[] element_list = line.split("\t");
				String sampleID = element_list[0];
				Float minDepth = Float.parseFloat(element_list[1]);
				Float maxDepth = Float.parseFloat(element_list[2]);
				ArrayList<Float> depthList = new ArrayList<Float>();
				depthList.add(minDepth);
				depthList.add(maxDepth);
				if (samples.contains(sampleID)){
					totalMinDepth += minDepth;
					totalMaxDepth += maxDepth;
					noOfSamplesWithDepth += 1;
				}
				if (element_list.length > 2) {
					Float minQUAL = Float.parseFloat(element_list[3]);
					depthList.add(minQUAL);
				} else {
					depthList.add((float) -1.0); // -1.0 as minGQ
				}
				sampleID2minMaxDepth.put(sampleID, depthList);
				lineNo++;
			}
			logger.info(noOfSamplesWithDepth + " samples from the input VCF.");
			logger.info("totalMinDepth= " + totalMinDepth + ".");
			logger.info("totalMaxDepth= " + totalMaxDepth + ".");
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Error: " + e.getMessage());
		}
		
		// 2012.5.1 make sure every sample from input VCF is in sampleMinMaxDepthFile.
		Iterator<String> samples_itr = samples.iterator();
		while (samples_itr.hasNext()) {
			String sampleID = (String) samples_itr.next();
			if (!sampleID2minMaxDepth.containsKey(sampleID)) {
				logger.error("sample " + sampleID + " from input VCF not available in sampleMinMaxDepthFile.\n");
				System.err.println("sample " + sampleID + " from input VCF not available in sampleMinMaxDepthFile.\n");
				System.exit(1);
			}
		}
		
		
		if (KEEP_ORIGINAL_CHR_COUNTS) {
			headerLines.add(new VCFFormatHeaderLine("AC_Orig", 1,
					VCFHeaderLineType.Integer, "Original AC"));
			headerLines.add(new VCFFormatHeaderLine("AF_Orig", 1,
					VCFHeaderLineType.Float, "Original AF"));
			headerLines.add(new VCFFormatHeaderLine("AN_Orig", 1,
					VCFHeaderLineType.Integer, "Original AN"));
		}

		if (vcfWriter != null) {
			// 2011-11-11 borrow the header from the input VCF
			vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
			// VCFUtils.getHeaderFields(this.getToolkit(),
			// SampleUtils.getUniqueSamplesFromRods(this.getToolkit()) returns
			// empty.
		}

		if (siteStatFname != null) {
			try {
				siteStatWriter = new FileWriter(siteStatFname);
				siteStatOutputHeader.add("chromosome");
				siteStatOutputHeader.add("pos");
				siteStatOutputHeader.add("missingRate");
				siteStatOutputHeader.add("AC");
				siteStatOutputHeader.add("AF");
				siteStatOutputHeader.add("hetFraction");
				siteStatOutputHeader.add("DP");
				siteStatOutputHeader.add("MQ");
				siteStatOutputHeader.add("HWE");
				siteStatOutputHeader.add("VDB");
				siteStatOutputHeader.add("strandBias");
				siteStatOutputHeader.add("baseQBias");
				siteStatOutputHeader.add("mapQBias");
				siteStatOutputHeader.add("tailDistanceBias");

				for (int i = 0; i < siteStatOutputHeader.size(); i++) {
					siteStatWriter.append(siteStatOutputHeader.get(i));
					if (i < siteStatOutputHeader.size() - 1)
						siteStatWriter.append("\t");
					else
						siteStatWriter.append('\n');
				}

			} catch (IOException e) {
				e.printStackTrace();
				System.err.println("Error: " + e.getMessage());
			}
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
		/*
		 * 2011-11-11 line from the example in the 20linesaver slides doesn't
		 * work here
		 */
		// Collection<VariantContext> vcs = tracker.getVariantContexts(ref,
		// "variant", null, context.getLocation(), true, true);

		if (vcs == null || vcs.size() == 0) {
			return 0;
		}

		for (VariantContext vc : vcs) {
			no_of_total_genotypes += vc.getCalledChrCount();
			if (vc.getAttributeAsInt("DP", 0)>=totalMinDepth && vc.getAttributeAsInt("DP", 0)<=totalMaxDepth){
				no_of_included_genotypes += vc.getCalledChrCount();
				
				/*
				GenotypesContext thisGC = vc.getGenotypes();
				// Map <String, Genotype> sampleID2genotype = vc.getGenotypes();
				ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
				Set<String> selectedSamples = new HashSet<String>(0);
				for (Genotype genotype : thisGC) {
					no_of_total_genotypes ++ ;
					String sampleID = genotype.getSampleName();
					if (sampleID2minMaxDepth.containsKey(sampleID)) {
						List<Float> depthList = sampleID2minMaxDepth.get(sampleID);
						Double genotypeGQ;
						Integer genotypeDepth;
						if (genotype.isCalled()) {
							genotypeGQ = genotype.getPhredScaledQual();
							// System.err.println(genotypeGQ.toString());
							// System.err.println(genotype.toString());
							genotypeDepth = getGenotypeDepth(genotype);
							// System.err.println(genotypeDepth.toString());
						} else {
							genotypeGQ = 0.0;
							genotypeDepth = 0;
						}
						if (genotypeDepth >= depthList.get(0)
								&& genotypeDepth <= depthList.get(1)
								&& genotypeGQ >= depthList.get(2)) {
							selectedSamples.add(sampleID);
							genotypes.add(genotype);
							no_of_included_genotypes++;
						}
					}
				}
				*/
				//VariantContext sub = vc.subContextFromSamples(selectedSamples);
				// VariantContext sub = vc.subContextFromGenotypes(genotypes,
				// vc.getAlleles());
				
				if (onlyKeepBiAllelicSNP == true) {
					// 2011.11.28
					if (vc.getAlleles().size() != 2) {
						no_of_nonBiAllelicLoci++;
						continue;
					}
				}
				if (vc.getGenotypes().size() > 0) { // output only when not all are
														// missing
					no_of_retained_loci++;
					//sub = subsetRecord(vc, selectedSamples);
					if (vcfWriter != null) {
						vcfWriter.add(vc);
					}
				}
			}
			if (siteStatWriter != null) {
				try {

					Vector<String> statVector = new Vector<String>();
					statVector.add(vc.getChr());

					Integer startPos = vc.getStart();
					statVector.add(startPos.toString());

					Float sndMissingRate = (vc.getGenotypes().size() - vc.getGenotypes().size()) / (float) vc.getNSamples();
					statVector.add(sndMissingRate.toString());

					/*
					 * 2011.12.20 sub.getNoCallCount() is 0.
					 * siteStatWriter.append("\t"); Float missingRate =
					 * sub.getNoCallCount()/(float)vc.getNSamples();
					 * siteStatWriter.write(missingRate.toString());
					 */
					statVector.add(getVariantAC(vc).toString());

					statVector.add(getVariantAF(vc).toString());

					Double hetFraction = vc.getHetCount()
							/ (double) vc.getNSamples();
					statVector.add(hetFraction.toString());

					Integer DP = vc.getAttributeAsInt("DP", -1);
					statVector.add(DP.toString());

					Integer mq = vc.getAttributeAsInt("MQ", -1);
					statVector.add(mq.toString());

					statVector.add(vc.getAttributeAsString("HWE", "-1"));

					statVector.add(vc.getAttributeAsString("VDB", "-1"));

					String pv4 = vc.getAttributeAsString("PV4", "");
					String[] pv4List = new String[] { "-1", "-1", "-1", "-1" };
					if (pv4.length() > 0) {
						pv4List = pv4.substring(1, pv4.length() - 1).split(",");
					}
					for (int i = 0; i < pv4List.length; i++) {
						statVector.add(pv4List[i]);
					}

					for (int i = 0; i < statVector.size(); i++) {
						siteStatWriter.write(statVector.get(i));
						if (i < statVector.size() - 1)
							siteStatWriter.append("\t");
						else
							siteStatWriter.append('\n');
					}

					/*
					 * siteStatWriter.append("\t"); Integer noOfGenotypes =
					 * vc.getGenotypes().size();
					 * siteStatWriter.write(noOfGenotypes.toString());
					 * 
					 * siteStatWriter.append("\t"); Integer noOfSamples =
					 * vc.getNSamples();
					 * siteStatWriter.write(noOfSamples.toString());
					 * 
					 * siteStatWriter.append("\t"); Integer
					 * noOfGenotypesAfterFilter = sub.getGenotypes().size();
					 * siteStatWriter
					 * .write(noOfGenotypesAfterFilter.toString());
					 */
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					System.err.println("Error: " + e.getMessage());
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
		logger.info(no_of_retained_loci + " loci retained.");
		logger.info("Number of total genotypes: " + no_of_total_genotypes);
		logger.info("Number of outputted genotypes: "
				+ no_of_included_genotypes);
		logger.info("Number of non-bi-allelic loci after GQ and Depth filter: "
				+ no_of_nonBiAllelicLoci);
		//2012. close the file, otherwise it might end up with some missing end (buffer not flushed to disk)
		if (siteStatWriter != null) {
			try {
				siteStatWriter.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.err.println("Error: " + e.getMessage());
			}
		}
	}

	/**
	 * 2011-11-11
	 * 	re-factored from some code in gatk/walkers/variantutils/SelectVariants.java
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

	/**
	 * 2011-12.20
	 * 
	 * @param genotype
	 * @return
	 */
	public Integer getVariantAC(VariantContext vc) {
		int AC = -1;
		/*
		 * if (vc.hasAttribute("AC") && vc.getAttribute("AC") instanceof
		 * Integer) { AC = vc.getAttributeAsInt("AC", 0); } else if
		 * (vc.hasAttribute("AC1") && vc.getAttribute("AC1") instanceof Integer)
		 * { AC = vc.getAttributeAsInt("AC1", 0); } else
		 */
		if (vc.isVariant()) {
			for (Allele allele : vc.getAlternateAlleles())
				AC = Math.max(AC, vc.getCalledChrCount(allele));
		} else
			// by default, the site is considered monomorphic
			AC = 0;

		return AC;
	}

	/**
	 * 2011-12.20
	 * 
	 * @param genotype
	 * @return
	 */
	public Double getVariantAF(VariantContext vc) {
		Double AF = -1.0;
		/*
		 * if (vc.hasAttribute("AF") && vc.getAttribute("AF") instanceof Double)
		 * { AF = vc.getAttributeAsDouble("AF", 0.0); } else if
		 * (vc.hasAttribute("AF1") && vc.getAttribute("AF1") instanceof Double)
		 * { AF = vc.getAttributeAsDouble("AF1", 0.0); } else
		 */
		if (vc.isVariant()) {
			for (Allele allele : vc.getAlternateAlleles()) {
				Double newAF = vc.getCalledChrCount(allele)
						/ (2 * (double) (vc.getNSamples() - vc.getNoCallCount()));
				AF = Math.max(AF, newAF);
			}
		} else
			// by default, the site is considered monomorphic
			AF = 0.0;

		return AF;
	}

	/**
	 * based on walker/variantutils/SelectVariants.java's subsetRecord() Helper
	 * method to modify some metadata stored in the INFO field (i.e. AN, AC,
	 * AF).
	 * 
	 * @param vc
	 *            the VariantContext record to subset
	 * @return the subsetted VariantContext
	 */
	public VariantContext subsetRecord(VariantContext vc, Set<String> samples) {
		if (samples == null || samples.isEmpty())
			return vc;

		final VariantContext sub = vc.subContextFromSamples(samples,
				vc.getAlleles());
		VariantContextBuilder builder = new VariantContextBuilder(sub);

		GenotypesContext newGC = sub.getGenotypes();

		// if we have fewer alternate alleles in the selected VC than in the
		// original VC, we need to strip out the GL/PLs (because they are no
		// longer accurate)
		if (vc.getAlleles().size() != sub.getAlleles().size())
			newGC = VariantContextUtils.stripPLs(sub.getGenotypes());

		builder.genotypes(newGC);
		int depth = 0;
		for (String sample : sub.getSampleNames()) {
			Genotype g = sub.getGenotype(sample);

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

		return builder.make();

	}

}
