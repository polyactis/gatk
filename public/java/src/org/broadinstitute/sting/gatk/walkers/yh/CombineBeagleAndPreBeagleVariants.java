/*
 * 2013.06.27 Copyright 2013, Yu Huang
 */
package org.broadinstitute.sting.gatk.walkers.yh;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.variantutils.CombineVariants;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.AlleleMapper;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.FilteredRecordMergeType;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.GenotypeMergeType;

/**
 * Combines VCF records from Beagle (v4) VCF and Pre-Beagle VCF files.
 * 	a child of GATK combineVariants to combine Beagle and Pre-Beagle VCF, taking priority of Beagle VCF, but add attributes from Pre-Beagle VCF to it.
 * 	It calculates PL (likelihood) and GQ (Genotype Quality) from Beagle GP (Genotype probability), instead of using old ones from Pre-Beagle VCF.
 *
 * <p>
 * 
 *
 * <h2>Input</h2>
 * <p>
 * Two VCFs with same set of samples and loci.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A combined VCF, with genotype and its attributes from 1st-priority file, and additional (non-existent in 1st-priority file) attributes from 2nd-priority VCF file.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * 
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CombineBeagleAndPreBeagleVariants \
 *   --variant:foo Beagle.vcf \
 *   --variant:bar input.vcf \
 *   -o output.vcf \
 *   -genotypeMergeOptions PRIORITIZE
 *   -priority foo,bar
 * </pre>
 *
 */

public class CombineBeagleAndPreBeagleVariants extends CombineVariants{

	
	private static void mergeGenotypes(GenotypesContext mergedGenotypes,
			VariantContext oneVC, AlleleMapper alleleMapping,
			boolean uniqifySamples) {
		for (Genotype g : oneVC.getGenotypes()) {
			String name = VariantContextUtils.mergedSampleName(oneVC.getSource(),
					g.getSampleName(), uniqifySamples);
			if (!mergedGenotypes.containsSample(name)) {
				// only add if the name is new
				Genotype newG = g;

				if (uniqifySamples || alleleMapping.needsRemapping()) {
					final List<Allele> alleles = alleleMapping.needsRemapping() ? alleleMapping
							.remap(g.getAlleles()) : g.getAlleles();
					newG = new Genotype(name, alleles, g.getLog10PError(),
							g.getFilters(), g.getAttributes(), g.isPhased());
				}

				mergedGenotypes.add(newG);
			}
		}
	}

	private static void addGenotypeAttributesFromOtherVariantContext(GenotypesContext mergedGenotypes,
			VariantContext oneVC, AlleleMapper alleleMapping) {
		/*
		 * 2013.09.15 deal with situation that newG (from oneVC, pre-Beagle genotype) could be "./." (instead of "./.:.:.") when it comes out
		 * 	of GATK and its depth is 0.
		 * 
		 */
		Double sndPriorityGenotypeLog10PError = null;
		double errorProbability = 0.0;
		double maxGenotypeProbability = -1;
		int indexOfGenotypeWithMaxProbability = -1;
		Vector<Double> genotypeProbabilityVector = new Vector<Double>();
		/*
		//2013.09.15 uniqueGenotypeAttributeNameSet is to store genotype FORMAT attribute names
		// so that for missing pre-Beagle genotypes ("./."), there is a way to add 
		Set<String> sndPriorityGenotypeAttributeKeySet = new HashSet<String>();
		boolean needLog10PError=false;
		boolean needPL = false;
		for (Genotype sndPriorityGenotype : oneVC.getGenotypes()) {
			sndPriorityGenotypeAttributeKeySet.addAll(sndPriorityGenotype.getAttributes().keySet());
			if (sndPriorityGenotype.hasLog10PError()){
				needLog10PError = true;
			}
			if (sndPriorityGenotype.hasLikelihoods()){
				needPL = true;
			}
		}
		logger.info("FORMAT at " + oneVC.getChr() + " : " + oneVC.getStart() + " is " + sndPriorityGenotypeAttributeKeySet.toString() + " \n");
		*/
		for (Genotype sndPriorityGenotype : oneVC.getGenotypes()) {
			String sampleName = sndPriorityGenotype.getSampleName();
			//sndPriorityGenotypeAttributeKeySet.addAll(sndPriorityGenotype.getAttributes().keySet());
			if (mergedGenotypes.containsSample(sampleName)) {
				Genotype oldG = mergedGenotypes.get(sampleName);
				//reset
				String[] genotypeProbabilityStrArray = oldG.getAttributeAsString("GP", "").split(",");
				//logger.info("GP at " + sndPriorityGenotype.getSampleName() + " is " + oldG.getAttributeAsString("GP", "") + " \n");
				double[] log10Likelihoods = new double[genotypeProbabilityStrArray.length];
				//double[] log10Likelihoods = new double[3];	//3 could cause Exception when there are >1 alternative alleles, etc.
				errorProbability = 0.0;
				indexOfGenotypeWithMaxProbability = -1;
				maxGenotypeProbability = -1;
				genotypeProbabilityVector.clear();
				
				// adding attributes to existing genotypes in mergedGenotypes
				if (alleleMapping.needsRemapping()) {
					final List<Allele> alleles = alleleMapping.needsRemapping() ? alleleMapping
							.remap(sndPriorityGenotype.getAlleles()) : sndPriorityGenotype.getAlleles();
					
					//sndPriorityGenotype = new Genotype(sampleName, alleles, g.getLog10PError(),
					//		g.getFilters(), g.getAttributes(), g.isPhased());
				}
				//get all attributes from sndPriorityGenotype.
				// 2013.09.15 
				Map<String, Object> newGenotypeAttributes = new TreeMap<String, Object>();
				//logger.info("PL at " + sndPriorityGenotype.getSampleName() + " is " + sndPriorityGenotype.getAttribute("PL") + " \n");
				// for missing genotype ("./."), PL is simply "." and sndPriorityGenotype.getAttributes() returns all attributes in FORMAT.
				for (final Map.Entry<String, Object> p: sndPriorityGenotype.getAttributes().entrySet()) {
					String attributeKey = p.getKey();
					if (attributeKey!="GT" && attributeKey!="PL" && attributeKey!="GQ"){
						//skip likelihoods
						//newGenotypeAttributes.put(attributeKey, ".,.,.");	//2103.09.15 put NA there if it does not have an attribute
						newGenotypeAttributes.put(attributeKey, sndPriorityGenotype.getAttribute(attributeKey, "."));
					}
				}
				// add/overwrite newGenotypeAttributes with attributes from oldG
				//newGenotypeAttributes.putAll(oldG.getAttributes());	//putAll results in java.lang.UnsupportedOperationException
				for (final Map.Entry<String, Object> p : oldG.getAttributes().entrySet()) {
					String key = p.getKey();
					newGenotypeAttributes.put(key, p.getValue());
					/*
					if (newGenotypeAttributes.containsKey(key)){
						newGenotypeAttributes.remove(key);
					}
					*/
				}
				//String[] genotypeProbabilityStrArray = oldG.getAttributeAsString("GP", "").split(",");
				// 2013.09.16 infer likelihood from genotype probability
				int genotypeIndex=0;
				for (final String genotypeProbabilityStr: genotypeProbabilityStrArray){
					Double genotypeProbability = Double.parseDouble(genotypeProbabilityStr);
					genotypeProbabilityVector.add(genotypeProbability);
					if (genotypeProbability>maxGenotypeProbability){
						maxGenotypeProbability = genotypeProbability;
						indexOfGenotypeWithMaxProbability  = genotypeIndex;
					}
					if (genotypeProbability<=0){
						log10Likelihoods[genotypeIndex] = -10.0;	//close to 0
					}
					else{
						log10Likelihoods[genotypeIndex] = Math.log10(genotypeProbability);
						
					}
					genotypeIndex++;
					
				}
				errorProbability = 1-maxGenotypeProbability;
				sndPriorityGenotypeLog10PError = Math.log10(errorProbability);	//could handle 0
				newGenotypeAttributes.remove("PL");	//remove PL , so that log10Likelihoods could overwrite it
				/*
				double oldGQ = oldG.getPhredScaledQual();
				if (!oldG.hasLog10PError()){
					//double GQ = sndPriorityGenotype.getLog10PError();
					double GQ = sndPriorityGenotype.getPhredScaledQual();
				}
				*/
				
				//now generate a new genotype using oldG's genotype and combined newGenotypeAttributes
				// note, log10PError could also be copied from sndPriorityGenotype if oldG doesn't have it
				// 2013.09.15 fake a sndPriorityGenotypeLog10PError if it's needed and sndPriorityGenotype (GATK missing genotype) not having it
				/*
				if (needLog10PError){
					if (sndPriorityGenotype.hasLog10PError()){
						sndPriorityGenotypeLog10PError = sndPriorityGenotype.getLog10PError();
					}
					else{
						//random default,to make GATK -> TrioCaller work
						sndPriorityGenotypeLog10PError = -5.0;
						//sndPriorityGenotypeLog10PError = null;	//would cause java.lang.NullPointerException
					}
				}
				*/
				oldG =  new Genotype(oldG.getSampleName(), oldG.getAlleles(), 
						oldG.hasLog10PError() ? oldG.getLog10PError() : sndPriorityGenotypeLog10PError, 
						oldG.filtersWereApplied() ? oldG.getFilters() : null,
						newGenotypeAttributes, oldG.isPhased(), log10Likelihoods);
				mergedGenotypes.replace(oldG);	//not sure if this does what i want
				//mergedGenotypes.add(sndPriorityGenotype);
			}
		}
	}
	
	/**
	 * Merges VariantContexts into a single hybrid. Takes genotypes for common
	 * samples in priority order, if provided. If uniquifySamples is true, the
	 * priority order is ignored and names are created by concatenating the VC
	 * name with the sample name
	 * 
	 * @param genomeLocParser
	 *            loc parser
	 * @param unsortedVCs
	 *            collection of unsorted VCs
	 * @param priorityListOfVCs
	 *            priority list detailing the order in which we should grab the
	 *            VCs
	 * @param filteredRecordMergeType
	 *            merge type for filtered records
	 * @param genotypeMergeOptions
	 *            merge option for genotypes
	 * @param annotateOrigin
	 *            should we annotate the set it came from?
	 * @param printMessages
	 *            should we print messages?
	 * @param setKey
	 *            the key name of the set
	 * @param filteredAreUncalled
	 *            are filtered records uncalled?
	 * @param mergeInfoWithMaxAC
	 *            should we merge in info from the VC with maximum allele count?
	 * @return new VariantContext representing the merge of unsortedVCs
	 */
	public static VariantContext simpleMerge(
			final GenomeLocParser genomeLocParser,
			final Collection<VariantContext> unsortedVCs,
			final List<String> priorityListOfVCs,
			final FilteredRecordMergeType filteredRecordMergeType,
			final GenotypeMergeType genotypeMergeOptions,
			final boolean annotateOrigin, final boolean printMessages,
			final String setKey, final boolean filteredAreUncalled,
			final boolean mergeInfoWithMaxAC) {
		if (unsortedVCs == null || unsortedVCs.size() == 0)
			return null;

		if (annotateOrigin && priorityListOfVCs == null)
			throw new IllegalArgumentException(
					"Cannot merge calls and annotate their origins without a complete priority list of VariantContexts");

		if (genotypeMergeOptions == GenotypeMergeType.REQUIRE_UNIQUE)
			VariantContextUtils.verifyUniqueSampleNames(unsortedVCs);

		final List<VariantContext> prepaddedVCs = VariantContextUtils.sortVariantContextsByPriority(
				unsortedVCs, priorityListOfVCs, genotypeMergeOptions);
		// Make sure all variant contexts are padded with reference base in case
		// of indels if necessary
		final List<VariantContext> VCs = new ArrayList<VariantContext>();

		for (final VariantContext vc : prepaddedVCs) {
			// also a reasonable place to remove filtered calls, if needed
			if (!filteredAreUncalled || vc.isNotFiltered())
				VCs.add(VariantContextUtils.createVariantContextWithPaddedAlleles(vc, false));
		}
		if (VCs.size() == 0) // everything is filtered out and we're
								// filteredAreUncalled
			return null;

		// establish the baseline info from the first VC
		final VariantContext first = VCs.get(0);
		final String name = first.getSource();
		final Allele refAllele = VariantContextUtils.determineReferenceAllele(VCs);
		Byte referenceBaseForIndel = null;

		final Set<Allele> alleles = new LinkedHashSet<Allele>();
		final Set<String> filters = new TreeSet<String>();
		final Map<String, Object> attributes = new TreeMap<String, Object>();
		final Set<String> inconsistentAttributes = new HashSet<String>();
		final Set<String> variantSources = new HashSet<String>(); // contains
																	// the set
																	// of
																	// sources
																	// we found
																	// in our
																	// set of
																	// VCs that
																	// are
																	// variant
		final Set<String> rsIDs = new LinkedHashSet<String>(1); // most of the
																// time there's
																// one id

		GenomeLoc loc = VariantContextUtils.getLocation(genomeLocParser, first);
		int depth = 0;
		int maxAC = -1;
		final Map<String, Object> attributesWithMaxAC = new TreeMap<String, Object>();
		double log10PError = 1;
		VariantContext vcWithMaxAC = null;
		GenotypesContext genotypes = GenotypesContext.create();

		// counting the number of filtered and variant VCs
		int nFiltered = 0;

		boolean remapped = false;

		// cycle through and add info from the other VCs, making sure the
		// loc/reference matches

		for (final VariantContext vc : VCs) {
			if (loc.getStart() != vc.getStart())
				throw new ReviewedStingException(
						"BUG: attempting to merge VariantContexts with different start sites: first="
								+ first.toString() + " second=" + vc.toString());

			if (VariantContextUtils.getLocation(genomeLocParser, vc).size() > loc.size())
				loc = VariantContextUtils.getLocation(genomeLocParser, vc); // get the longest
														// location

			nFiltered += vc.isFiltered() ? 1 : 0;
			if (vc.isVariant())
				variantSources.add(vc.getSource());

			AlleleMapper alleleMapping = VariantContextUtils.resolveIncompatibleAlleles(refAllele,
					vc, alleles);
			remapped = remapped || alleleMapping.needsRemapping();

			alleles.addAll(alleleMapping.values());

			mergeGenotypes(genotypes, vc, alleleMapping,
					genotypeMergeOptions == GenotypeMergeType.UNIQUIFY);

			log10PError = Math.min(log10PError,
					vc.isVariant() ? vc.getLog10PError() : 1);

			filters.addAll(vc.getFilters());

			if (referenceBaseForIndel == null)
				referenceBaseForIndel = vc.getReferenceBaseForIndel();

			//
			// add attributes
			//
			// special case DP (add it up) and ID (just preserve it)
			//
			if (vc.hasAttribute(VCFConstants.DEPTH_KEY))
				depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
			if (vc.hasID())
				rsIDs.add(vc.getID());
			if (mergeInfoWithMaxAC
					&& vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
				String rawAlleleCounts = vc.getAttributeAsString(
						VCFConstants.ALLELE_COUNT_KEY, null);
				// lets see if the string contains a , separator
				if (rawAlleleCounts
						.contains(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) {
					List<String> alleleCountArray = Arrays
							.asList(rawAlleleCounts.substring(1,
									rawAlleleCounts.length() - 1).split(
									VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
					for (String alleleCount : alleleCountArray) {
						final int ac = Integer.valueOf(alleleCount.trim());
						if (ac > maxAC) {
							maxAC = ac;
							vcWithMaxAC = vc;
						}
					}
				} else {
					final int ac = Integer.valueOf(rawAlleleCounts);
					if (ac > maxAC) {
						maxAC = ac;
						vcWithMaxAC = vc;
					}
				}
			}

			for (final Map.Entry<String, Object> p : vc.getAttributes()
					.entrySet()) {
				String key = p.getKey();
				// if we don't like the key already, don't go anywhere
				if (!inconsistentAttributes.contains(key)) {
					final boolean alreadyFound = attributes.containsKey(key);
					final Object boundValue = attributes.get(key);
					final boolean boundIsMissingValue = alreadyFound
							&& boundValue.equals(VCFConstants.MISSING_VALUE_v4);

					if (alreadyFound && !boundValue.equals(p.getValue())
							&& !boundIsMissingValue) {
						// we found the value but we're inconsistent, put it in
						// the exclude list
						// System.out.printf("Inconsistent INFO values: %s => %s and %s%n",
						// key, boundValue, p.getValue());
						inconsistentAttributes.add(key);
						attributes.remove(key);
					} else if (!alreadyFound || boundIsMissingValue) { // no
																		// value
						// if ( vc != first )
						// System.out.printf("Adding key %s => %s%n",
						// p.getKey(), p.getValue());
						attributes.put(key, p.getValue());
					}
				}
			}
		}
		// 2013.06.27 Yu Huang	add genotype attributes from 2nd-priority VCF to 1st-priority genotype
		VariantContext secondPriorityVC = VCs.get(1);
		AlleleMapper alleleMapping = VariantContextUtils.resolveIncompatibleAlleles(refAllele,
				secondPriorityVC, alleles);
		addGenotypeAttributesFromOtherVariantContext(genotypes, secondPriorityVC, alleleMapping);

		// if we have more alternate alleles in the merged VC than in one or
		// more of the
		// original VCs, we need to strip out the GL/PLs (because they are no
		// longer accurate), as well as allele-dependent attributes like AC,AF
		for (final VariantContext vc : VCs) {
			if (vc.alleles.size() == 1)
				continue;
			if (VariantContextUtils.hasPLIncompatibleAlleles(alleles, vc.alleles)) {
				if (!genotypes.isEmpty())
					logger.debug(String
							.format("Stripping PLs at %s due incompatible alleles merged=%s vs. single=%s",
									genomeLocParser.createGenomeLoc(vc),
									alleles, vc.alleles));
				genotypes = VariantContextUtils.stripPLs(genotypes);
				// this will remove stale AC,AF attributed from vc
				VariantContextUtils.calculateChromosomeCounts(vc, attributes, true);
				break;
			}
		}

		// take the VC with the maxAC and pull the attributes into a modifiable
		// map
		if (mergeInfoWithMaxAC && vcWithMaxAC != null) {
			attributesWithMaxAC.putAll(vcWithMaxAC.getAttributes());
		}

		// if at least one record was unfiltered and we want a union, clear all
		// of the filters
		if ((filteredRecordMergeType == FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED && nFiltered != VCs
				.size())
				|| filteredRecordMergeType == FilteredRecordMergeType.KEEP_UNCONDITIONAL)
			filters.clear();

		if (annotateOrigin) { // we care about where the call came from
			String setValue;
			if (nFiltered == 0
					&& variantSources.size() == priorityListOfVCs.size()) // nothing
																			// was
																			// unfiltered
				setValue = VariantContextUtils.MERGE_INTERSECTION;
			else if (nFiltered == VCs.size()) // everything was filtered out
				setValue = VariantContextUtils.MERGE_FILTER_IN_ALL;
			else if (variantSources.isEmpty()) // everyone was reference
				setValue = VariantContextUtils.MERGE_REF_IN_ALL;
			else {
				final LinkedHashSet<String> s = new LinkedHashSet<String>();
				for (final VariantContext vc : VCs)
					if (vc.isVariant())
						s.add(vc.isFiltered() ? VariantContextUtils.MERGE_FILTER_PREFIX
								+ vc.getSource() : vc.getSource());
				setValue = Utils.join("-", s);
			}

			if (setKey != null) {
				attributes.put(setKey, setValue);
				if (mergeInfoWithMaxAC && vcWithMaxAC != null) {
					attributesWithMaxAC.put(setKey, vcWithMaxAC.getSource());
				}
			}
		}

		if (depth > 0)
			attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));

		final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils
				.join(",", rsIDs);

		final VariantContextBuilder builder = new VariantContextBuilder()
				.source(name).id(ID);
		builder.loc(loc.getContig(), loc.getStart(), loc.getStop());
		builder.alleles(alleles);
		builder.genotypes(genotypes);
		builder.log10PError(log10PError);
		builder.filters(filters).attributes(
				mergeInfoWithMaxAC ? attributesWithMaxAC : attributes);
		builder.referenceBaseForIndel(referenceBaseForIndel);

		// Trim the padded bases of all alleles if necessary
		final VariantContext merged = VariantContextUtils.createVariantContextWithTrimmedAlleles(builder
				.make());
		if (printMessages && remapped)
			System.out.printf("Remapped => %s%n", merged);
		return merged;
	}

	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		if (tracker == null) // RodWalkers can make funky map calls
			return 0;

		// get all of the vcf rods at this locus
		// Need to provide reference bases to simpleMerge starting at current
		// locus
		Collection<VariantContext> vcs = tracker.getValues(variants,
				context.getLocation());

		if (sitesOnlyVCF) {
			vcs = VariantContextUtils.sitesOnlyVariantContexts(vcs);
		}

		if (ASSUME_IDENTICAL_SAMPLES) {
			for (final VariantContext vc : vcs) {
				vcfWriter.add(vc);
			}

			return vcs.isEmpty() ? 0 : 1;
		}

		int numFilteredRecords = 0;
		for (VariantContext vc : vcs) {
			if (vc.filtersWereApplied() && vc.isFiltered())
				numFilteredRecords++;
		}

		if (minimumN > 1 && (vcs.size() - numFilteredRecords < minimumN))
			return 0;

		List<VariantContext> mergedVCs = new ArrayList<VariantContext>();

		if (multipleAllelesMergeType == VariantContextUtils.MultipleAllelesMergeType.BY_TYPE) {
			Map<VariantContext.Type, List<VariantContext>> VCsByType = VariantContextUtils
					.separateVariantContextsByType(vcs);

			// TODO -- clean this up in a refactoring
			// merge NO_VARIATION into another type of variant (based on the
			// ordering in VariantContext.Type)
			if (VCsByType.containsKey(VariantContext.Type.NO_VARIATION)
					&& VCsByType.size() > 1) {
				final List<VariantContext> refs = VCsByType
						.remove(VariantContext.Type.NO_VARIATION);
				for (VariantContext.Type type : VariantContext.Type.values()) {
					if (VCsByType.containsKey(type)) {
						VCsByType.get(type).addAll(refs);
						break;
					}
				}
			}

			// iterate over the types so that it's deterministic
			for (VariantContext.Type type : VariantContext.Type.values()) {
				if (VCsByType.containsKey(type))
					mergedVCs.add(simpleMerge(getToolkit()
							.getGenomeLocParser(), VCsByType.get(type),
							priority, filteredRecordsMergeType,
							genotypeMergeOption, true, printComplexMerges,
							SET_KEY, filteredAreUncalled,
							MERGE_INFO_WITH_MAX_AC));
			}
		} else if (multipleAllelesMergeType == VariantContextUtils.MultipleAllelesMergeType.MIX_TYPES) {
			mergedVCs.add(simpleMerge(getToolkit()
					.getGenomeLocParser(), vcs, priority,
					filteredRecordsMergeType, genotypeMergeOption, true,
					printComplexMerges, SET_KEY, filteredAreUncalled,
					MERGE_INFO_WITH_MAX_AC));
		} else {
			logger.warn("Ignoring all records at site " + ref.getLocus());
		}

		for (VariantContext mergedVC : mergedVCs) {
			// only operate at the start of events
			if (mergedVC == null)
				continue;

			final VariantContextBuilder builder = new VariantContextBuilder(
					mergedVC);
			// re-compute chromosome counts
			VariantContextUtils.calculateChromosomeCounts(builder, false);
			if (minimalVCF)
				VariantContextUtils.pruneVariantContext(builder,
						Arrays.asList(SET_KEY));
			vcfWriter.add(builder.make());
		}

		return vcs.isEmpty() ? 0 : 1;
	}
}
