/*
 * Copyright (c) 2010.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class IndelGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {
    private final int HAPLOTYPE_SIZE;

    private final int minIndelCountForGenotyping;
    private final boolean getAlleleListFromVCF;

    private boolean DEBUG = false;
    private final boolean doMultiAllelicCalls = true;
    private boolean ignoreSNPAllelesWhenGenotypingIndels = false;
    private PairHMMIndelErrorModel pairModel;

    private static ThreadLocal<HashMap<PileupElement, LinkedHashMap<Allele, Double>>> indelLikelihoodMap =
            new ThreadLocal<HashMap<PileupElement, LinkedHashMap<Allele, Double>>>() {
                protected synchronized HashMap<PileupElement, LinkedHashMap<Allele, Double>> initialValue() {
                    return new HashMap<PileupElement, LinkedHashMap<Allele, Double>>();
                }
            };

    private LinkedHashMap<Allele, Haplotype> haplotypeMap;

    // gdebug removeme
    // todo -cleanup
    private GenomeLoc lastSiteVisited;
    private ArrayList<Allele> alleleList;

    static {
        indelLikelihoodMap.set(new HashMap<PileupElement, LinkedHashMap<Allele, Double>>());
    }


    public IndelGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC, logger);
        pairModel = new PairHMMIndelErrorModel(UAC.INDEL_GAP_OPEN_PENALTY, UAC.INDEL_GAP_CONTINUATION_PENALTY,
                UAC.OUTPUT_DEBUG_INDEL_INFO, !UAC.DONT_DO_BANDED_INDEL_COMPUTATION);
        alleleList = new ArrayList<Allele>();
        getAlleleListFromVCF = UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
        minIndelCountForGenotyping = UAC.MIN_INDEL_COUNT_FOR_GENOTYPING;
        HAPLOTYPE_SIZE = UAC.INDEL_HAPLOTYPE_SIZE;
        DEBUG = UAC.OUTPUT_DEBUG_INDEL_INFO;

        haplotypeMap = new LinkedHashMap<Allele, Haplotype>();
        ignoreSNPAllelesWhenGenotypingIndels = UAC.IGNORE_SNP_ALLELES;
    }


    public static ArrayList<Allele> computeConsensusAlleles(ReferenceContext ref,
                                                      Map<String, AlignmentContext> contexts,
                                                      AlignmentContextUtils.ReadOrientation contextType, GenomeLocParser locParser,
                                                      int minIndelCountForGenotyping, boolean doMultiAllelicCalls) {
        Allele refAllele = null, altAllele = null;
        GenomeLoc loc = ref.getLocus();
        ArrayList<Allele> aList = new ArrayList<Allele>();

        HashMap<String, Integer> consensusIndelStrings = new HashMap<String, Integer>();

        int insCount = 0, delCount = 0;
        // quick check of total number of indels in pileup
        for (Map.Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadBackedExtendedEventPileup indelPileup = context.getExtendedEventPileup();
            insCount += indelPileup.getNumberOfInsertions();
            delCount += indelPileup.getNumberOfDeletions();
        }

        if (insCount < minIndelCountForGenotyping && delCount < minIndelCountForGenotyping)
            return aList;

        for (Map.Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            // todo -- warning, can be duplicating expensive partition here
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadBackedExtendedEventPileup indelPileup = context.getExtendedEventPileup();


            for (ExtendedEventPileupElement p : indelPileup.toExtendedIterable()) {
                //SAMRecord read = p.getRead();
                GATKSAMRecord read = ReadClipper.hardClipAdaptorSequence(p.getRead());
                if (read == null)
                    continue;
                if (ReadUtils.is454Read(read)) {
                    continue;
                }

/*                if (DEBUG && p.isIndel()) {
                    System.out.format("Read: %s, cigar: %s, aln start: %d, aln end: %d, p.len:%d, Type:%s, EventBases:%s\n",
                            read.getReadName(),read.getCigar().toString(),read.getAlignmentStart(),read.getAlignmentEnd(),
                            p.getEventLength(),p.getType().toString(), p.getEventBases());
                }
   */

                String indelString = p.getEventBases();
                if (p.isInsertion()) {
                    boolean foundKey = false;
                    // copy of hashmap into temp arrayList
                    ArrayList<Pair<String,Integer>> cList = new ArrayList<Pair<String,Integer>>();
                    for (String s : consensusIndelStrings.keySet()) {
                        cList.add(new Pair<String, Integer>(s,consensusIndelStrings.get(s)));
                    }

                    if (read.getAlignmentEnd() == loc.getStart()) {
                        // first corner condition: a read has an insertion at the end, and we're right at the insertion.
                        // In this case, the read could have any of the inserted bases and we need to build a consensus

                        for (int k=0; k < cList.size(); k++) {
                            String s = cList.get(k).getFirst();
                            int cnt = cList.get(k).getSecond();
                            // case 1: current insertion is prefix of indel in hash map
                            if (s.startsWith(indelString)) {
                                cList.set(k,new Pair<String, Integer>(s,cnt+1));
                                foundKey = true;
                            }
                            else if (indelString.startsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                foundKey = true;
                                cList.set(k,new Pair<String, Integer>(indelString,cnt+1));
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                            cList.add(new Pair<String, Integer>(indelString,1));

                    }
                    else if (read.getAlignmentStart() == loc.getStart()+1) {
                        // opposite corner condition: read will start at current locus with an insertion
                        for (int k=0; k < cList.size(); k++) {
                            String s = cList.get(k).getFirst();
                            int cnt = cList.get(k).getSecond();
                            if (s.endsWith(indelString)) {
                                // case 1: current insertion (indelString) is suffix of indel in hash map (s)
                                cList.set(k,new Pair<String, Integer>(s,cnt+1));
                                foundKey = true;
                            }
                            else if (indelString.endsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                foundKey = true;
                                cList.set(k,new Pair<String, Integer>(indelString,cnt+1));
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                            cList.add(new Pair<String, Integer>(indelString,1));


                    }
                    else {
                        // normal case: insertion somewhere in the middle of a read: add count to arrayList
                        int cnt = consensusIndelStrings.containsKey(indelString)? consensusIndelStrings.get(indelString):0;
                        cList.add(new Pair<String, Integer>(indelString,cnt+1));
                    }

                    // copy back arrayList into hashMap
                    consensusIndelStrings.clear();
                    for (Pair<String,Integer> pair : cList) {
                        consensusIndelStrings.put(pair.getFirst(),pair.getSecond());
                    }

                }
                else if (p.isDeletion()) {
                    indelString = String.format("D%d",p.getEventLength());
                    int cnt = consensusIndelStrings.containsKey(indelString)? consensusIndelStrings.get(indelString):0;
                    consensusIndelStrings.put(indelString,cnt+1);

                }
            }

        }

        Collection<VariantContext> vcs = new ArrayList<VariantContext>();
        int maxAlleleCnt = 0;
        String bestAltAllele = "";

        for (String s : consensusIndelStrings.keySet()) {
            int curCnt = consensusIndelStrings.get(s), stop = 0;
            // if observed count if above minimum threshold, we will genotype this allele
            if (curCnt < minIndelCountForGenotyping)
                continue;

            if (s.startsWith("D")) {
                // get deletion length
                int dLen = Integer.valueOf(s.substring(1));
                // get ref bases of accurate deletion
                int startIdxInReference = 1 + loc.getStart() - ref.getWindow().getStart();
                stop = loc.getStart() + dLen;
                byte[] refBases = Arrays.copyOfRange(ref.getBases(), startIdxInReference, startIdxInReference + dLen);

                if (Allele.acceptableAlleleBases(refBases)) {
                    refAllele = Allele.create(refBases, true);
                    altAllele = Allele.create(Allele.NULL_ALLELE_STRING, false);
                }
                else continue; // don't go on with this allele if refBases are non-standard
            } else {
                // insertion case
                if (Allele.acceptableAlleleBases(s)) {
                    refAllele = Allele.create(Allele.NULL_ALLELE_STRING, true);
                    altAllele = Allele.create(s, false);
                    stop = loc.getStart();
                }
                else continue; // go on to next allele if consensus insertion has any non-standard base.
            }


            ArrayList vcAlleles = new ArrayList<Allele>();
            vcAlleles.add(refAllele);
            vcAlleles.add(altAllele);

            final VariantContextBuilder builder = new VariantContextBuilder().source("");
            builder.loc(loc.getContig(), loc.getStart(), stop);
            builder.alleles(vcAlleles);
            builder.referenceBaseForIndel(ref.getBase());
            builder.noGenotypes();
            if (doMultiAllelicCalls)
                vcs.add(builder.make());
            else {
                if (curCnt > maxAlleleCnt) {
                    maxAlleleCnt = curCnt;
                    vcs.clear();
                    vcs.add(builder.make());
                }

            }
        }

        if (vcs.isEmpty())
            return aList; // nothing else to do, no alleles passed minimum count criterion

        VariantContext mergedVC = VariantContextUtils.simpleMerge(locParser, vcs, null, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, VariantContextUtils.GenotypeMergeType.UNSORTED, false, false, null, false, false);

        aList = new ArrayList<Allele>(mergedVC.getAlleles());

        return aList;

    }

    private final static EnumSet<VariantContext.Type> allowableTypes = EnumSet.of(VariantContext.Type.INDEL, VariantContext.Type.MIXED);

    public VariantContext getLikelihoods(final RefMetaDataTracker tracker,
                                         final ReferenceContext ref,
                                         final Map<String, AlignmentContext> contexts,
                                         final AlignmentContextUtils.ReadOrientation contextType,
                                         final GenotypePriors priors,
                                         final List<Allele> alternateAllelesToUse,
                                         final boolean useBAQedPileup,
                                         final GenomeLocParser locParser) {

        if (tracker == null)
            return null;

        GenomeLoc loc = ref.getLocus();
        Allele refAllele, altAllele;
        VariantContext vc = null;

        if (!ref.getLocus().equals(lastSiteVisited)) {
            // starting a new site: clear allele list
            alleleList.clear();
            lastSiteVisited = ref.getLocus();
            indelLikelihoodMap.set(new HashMap<PileupElement, LinkedHashMap<Allele, Double>>());
            haplotypeMap.clear();

            if (getAlleleListFromVCF) {
                for (final VariantContext vc_input : tracker.getValues(UAC.alleles, loc)) {
                    if (vc_input != null &&
                            allowableTypes.contains(vc_input.getType()) &&
                            ref.getLocus().getStart() == vc_input.getStart()) {
                        vc = vc_input;
                        break;
                    }
                }
                // ignore places where we don't have a variant
                if (vc == null)
                    return null;

                alleleList.clear();
                if (ignoreSNPAllelesWhenGenotypingIndels) {
                    // if there's an allele that has same length as the reference (i.e. a SNP or MNP), ignore it and don't genotype it
                    for (Allele a : vc.getAlleles())
                        if (a.isNonReference() && a.getBases().length == vc.getReference().getBases().length)
                            continue;
                        else
                            alleleList.add(a);

                } else {
                    for (Allele a : vc.getAlleles())
                        alleleList.add(a);
                }

            } else {
                alleleList = computeConsensusAlleles(ref, contexts, contextType, locParser, minIndelCountForGenotyping,doMultiAllelicCalls);
                if (alleleList.isEmpty())
                    return null;
            }
        }
        // protect against having an indel too close to the edge of a contig
        if (loc.getStart() <= HAPLOTYPE_SIZE)
            return null;

        // check if there is enough reference window to create haplotypes (can be an issue at end of contigs)
        if (ref.getWindow().getStop() < loc.getStop() + HAPLOTYPE_SIZE)
            return null;
        if (!(priors instanceof DiploidIndelGenotypePriors))
            throw new StingException("Only diploid-based Indel priors are supported in the INDEL GL model");

        if (alleleList.isEmpty())
            return null;

        refAllele = alleleList.get(0);
        altAllele = alleleList.get(1);

        // look for alt allele that has biggest length distance to ref allele
        int maxLenDiff = 0;
        for (Allele a : alleleList) {
            if (a.isNonReference()) {
                int lenDiff = Math.abs(a.getBaseString().length() - refAllele.getBaseString().length());
                if (lenDiff > maxLenDiff) {
                    maxLenDiff = lenDiff;
                    altAllele = a;
                }
            }
        }

        final int eventLength = altAllele.getBaseString().length() - refAllele.getBaseString().length();
        final int hsize = (int) ref.getWindow().size() - Math.abs(eventLength) - 1;
        final int numPrefBases = ref.getLocus().getStart() - ref.getWindow().getStart() + 1;

        if (hsize <= 0) {
            logger.warn(String.format("Warning: event at location %s can't be genotyped, skipping", loc.toString()));
            return null;
        }
        haplotypeMap = Haplotype.makeHaplotypeListFromAlleles(alleleList, loc.getStart(),
                ref, hsize, numPrefBases);

        // start making the VariantContext
        final int endLoc = calculateEndPos(alleleList, refAllele, loc);
        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), endLoc, alleleList).referenceBaseForIndel(ref.getBase());

        // create the genotypes; no-call everyone for now
        GenotypesContext genotypes = GenotypesContext.create();
        final List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        // For each sample, get genotype likelihoods based on pileup
        // compute prior likelihoods on haplotypes, and initialize haplotype likelihood matrix with them.

        for (Map.Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            ReadBackedPileup pileup = null;
            if (context.hasExtendedEventPileup())
                pileup = context.getExtendedEventPileup();
            else if (context.hasBasePileup())
                pileup = context.getBasePileup();

            if (pileup != null) {
                final double[] genotypeLikelihoods = pairModel.computeReadHaplotypeLikelihoods(pileup, haplotypeMap, ref, eventLength, getIndelLikelihoodMap());
                GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(genotypeLikelihoods);

                HashMap<String, Object> attributes = new HashMap<String, Object>();
                attributes.put(VCFConstants.DEPTH_KEY, getFilteredDepth(pileup));
                attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, likelihoods);
                genotypes.add(new Genotype(sample.getKey(), noCall, Genotype.NO_LOG10_PERROR, null, attributes, false));

                if (DEBUG) {
                    System.out.format("Sample:%s Alleles:%s GL:", sample.getKey(), alleleList.toString());
                    for (int k = 0; k < genotypeLikelihoods.length; k++)
                        System.out.format("%1.4f ", genotypeLikelihoods[k]);
                    System.out.println();
                }
            }
        }

        return builder.genotypes(genotypes).make();
    }

    private int calculateEndPos(Collection<Allele> alleles, Allele refAllele, GenomeLoc loc) {
        // for indels, stop location is one more than ref allele length
        boolean hasNullAltAllele = false;
        for (Allele a : alleles) {
            if (a.isNull()) {
                hasNullAltAllele = true;
                break;
            }
        }

        int endLoc = loc.getStart() + refAllele.length();
        if (!hasNullAltAllele)
            endLoc--;

        return endLoc;
    }

    public static HashMap<PileupElement, LinkedHashMap<Allele, Double>> getIndelLikelihoodMap() {
        return indelLikelihoodMap.get();
    }

    // Overload function in GenotypeLikelihoodsCalculationModel so that, for an indel case, we consider a deletion as part of the pileup,
    // so that per-sample DP will include deletions covering the event.
    protected int getFilteredDepth(ReadBackedPileup pileup) {
        int count = 0;
        for (PileupElement p : pileup) {
            if (p.isDeletion() || p.isInsertionAtBeginningOfRead() || BaseUtils.isRegularBase(p.getBase()))
                count++;
        }

        return count;
    }

}