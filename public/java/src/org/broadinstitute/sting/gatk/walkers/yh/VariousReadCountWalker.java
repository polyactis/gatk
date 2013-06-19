package org.broadinstitute.sting.gatk.walkers.yh;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/*
 * 2011-11-3 yh
 * 1. count the number of total reads aligned, number of total, percentage of aligned
 * 2. count the number of pairs mapped, the number of singletons (its mate un-mapped) mapped
 * 3. insert size statistic (mean/median/rmsd)
 */
@Requires({ DataSource.READS, DataSource.REFERENCE })
public class VariousReadCountWalker extends ReadWalker<Integer, Integer>  {

	@Output
	private PrintStream out;

	@Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads. Defaults to -1. " +
			"Filtering would result in the total number of reads smaller than that inferred from other numbers.",
			required = false)
	int minMappingQuality = -1;
	@Argument(fullName = "maxMappingQuality", doc = "Maximum mapping quality of reads. Defaults to 2^31-1 (Integer.MAX_VALUE).", 
			required = false)
	int maxMappingQuality = Integer.MAX_VALUE;
	/**
	 * How many reads are the first in a pair, based on flag 0x0040 from the SAM
	 * spec.
	 */
	private long firstOfPair = 0;

	/**
	 * How many reads are the second in a pair, based on flag 0x0080 from the
	 * SAM spec.
	 */
	private long secondOfPair = 0;

	private Long numberOfReadsAligned=(long) 0;
	private Long numberOfReads=(long) 0;
	private Long numberOfDuplicatedPairs=(long) 0;
	private Long numberOfOtherPairs = (long) 0;
	private Double numberOfPairsOnSameContig=0.0;
	private Double numberOfPairsOnDifferentContigs=0.0;
	private Long numberOfSingletonsMapped=(long) 0;
	
	private String separator = "\t";	//used in output to separate columns
	
	//Pair<GenomeLoc,DepthOfCoverageStats> targetStats = new Pair<GenomeLoc,DepthOfCoverageStats>(
	//		targetAggregator.first, targetAggregator.second.getCoverageByAggregationType(type));
	
	Map<String, Integer> readPairName2Count = new HashMap<String, Integer>();
	ArrayList<Integer> insertSizeLs = new ArrayList<Integer>(0);
	
	private Long cumulativeInsertSize = (long) 0;
	
	private String readGroup = "";
	private String referenceName = "";
	private Long referenceLength = (long) -1;
	
	// Do we actually want to operate on the context?
	/**
	 * Must return true for reads that need to be processed. Reads, for which
	 * this method return false will be skipped by the engine and never passed
	 * to the walker.
	 */
	@Override
	public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
		if (read.getMappingQuality()>=minMappingQuality){
			return true;
		}
		else{
			return false;
		}
	}

	public boolean readWellMapped(SAMRecord read) {
		/*
		 * check read.getMateUnmappedFlag() and read.getMappingQuality()>=minMappingQuality
		 */
		if (read.getReadPairedFlag() && !read.getMateUnmappedFlag() && read.getMappingQuality()>=minMappingQuality){
			return true;
		}
		else{
			return false;
		}
	}
	
	@Override
	public Integer map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {
		if (!read.getReadUnmappedFlag()){
			numberOfReadsAligned ++;
		}
		numberOfReads ++;
		
		if (readGroup==""){
			readGroup = read.getReadGroup().getReadGroupId();
		}
		
		if (read.getDuplicateReadFlag())
		{
			numberOfDuplicatedPairs += 1;
		}
		else{
			if (readWellMapped(read) && !read.getReadUnmappedFlag())	//both mates are mapped.
			{
				//String readPairName = read.getReadName();	//both mates have the same name. the trailing #0/1 or #0/2 is removed during alignment or whatever.
				/*	//out of memory on 4G limit when running on whole 30X bam
				if (readPairName2Count.containsKey(readPairName))
				{
					readPairName2Count.put(readPairName, readPairName2Count.get(readPairName)+1);
				}
				else{
					readPairName2Count.put(readPairName, 1);
				}
				int readPairCount = readPairName2Count.get(readPairName);
				*/
				//if (readPairCount==1){//count only once for every pair
					boolean isProperlyPaired = read.getProperPairFlag();
					//Note: the mate's mapq is not checked and some of them with low-mapq would be included as well.
					if (read.getReferenceIndex()==read.getMateReferenceIndex()){
						numberOfPairsOnSameContig = numberOfPairsOnSameContig + 1;	//0.5 because each pair is counted twice.
						//insertSizeLs.add(read.getInferredInsertSize());	//out of memory on 4G limit when running on whole 30X bam
						cumulativeInsertSize = cumulativeInsertSize + Math.abs((long) read.getInferredInsertSize());
					}
					else{
						numberOfPairsOnDifferentContigs = numberOfPairsOnDifferentContigs + 1;
					}
				//}
			}
			else if (!read.getReadPairedFlag() || (read.getMateUnmappedFlag() && !read.getReadUnmappedFlag())){	//its mate is not mapped.
				numberOfSingletonsMapped ++;
			}
			else{
				numberOfOtherPairs ++;
			}
		}
		if (referenceLength==-1){
			referenceLength = (long) read.getHeader().getSequence(read.getReferenceIndex()).getSequenceLength();
		}
		if (referenceName==""){
			referenceName = read.getReferenceName();
		}
		return 1;
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
		return value + sum;
	}
	@Override
	public Integer reduceInit() {
		return 0;
	}

	

	@Override
	public void onTraversalDone(Integer sum) {
		// TODO Auto-generated method stub
		// TODO generate the RMSD of inferred insert size list and its median (if easy)
		// TODO report the number properly paired
		StringBuilder summaryHeader = new StringBuilder();
		summaryHeader.append("readGroup");
		summaryHeader.append(separator);
		summaryHeader.append("firstReferenceName");
		summaryHeader.append(separator);
		summaryHeader.append("firstReferenceLength");
		summaryHeader.append(separator);
		summaryHeader.append("numberOfReads");
		summaryHeader.append(separator);
		summaryHeader.append("numberOfReadsAligned");
		summaryHeader.append(separator);
		summaryHeader.append("numberOfDuplicatedPairs");
		summaryHeader.append(separator);
		summaryHeader.append("numberOfSingletonsMapped");
		summaryHeader.append(separator);
		summaryHeader.append("numberOfPairsOnSameContig");
		summaryHeader.append(separator);
		summaryHeader.append("numberOfPairsOnDifferentContigs");
		summaryHeader.append(separator);
		summaryHeader.append("numberOfOtherPairs");
		summaryHeader.append(separator);
		summaryHeader.append("meanInsertSize");
		out.printf("%s%n", summaryHeader);
		
		StringBuilder targetSummary = new StringBuilder();
		targetSummary.append(readGroup);
		targetSummary.append(separator);
		targetSummary.append(referenceName);
		targetSummary.append(separator);
		targetSummary.append(referenceLength);
		targetSummary.append(separator);
		targetSummary.append(numberOfReads.toString());
		targetSummary.append(separator);
		targetSummary.append(numberOfReadsAligned.toString());
		targetSummary.append(separator);
		targetSummary.append(numberOfDuplicatedPairs.toString());
		targetSummary.append(separator);
		targetSummary.append(numberOfSingletonsMapped.toString());
		targetSummary.append(separator);
		targetSummary.append(String.format("%.1f", numberOfPairsOnSameContig));
		targetSummary.append(separator);
		targetSummary.append(String.format("%.1f", numberOfPairsOnDifferentContigs));
		targetSummary.append(separator);
		targetSummary.append(numberOfOtherPairs.toString());
		targetSummary.append(separator);
		if (numberOfPairsOnSameContig>0)
			targetSummary.append(String.format("%.1f", cumulativeInsertSize.floatValue()/numberOfPairsOnSameContig));
		else
			targetSummary.append(String.format("-1"));
		
		out.printf("%s%n", targetSummary);
		
		logger.info("[REDUCE RESULT] Traversal result is: " + sum);
	}

	/**
	 * General interval reduce routine called after all of the traversals are
	 * done
	 * 
	 * @param results
	 *            interval reduce results

	@Override
	public void onTraversalDone(List<Pair<GenomeLoc, Integer>> results) {
		for (Pair<GenomeLoc, ReduceType> result : results) {
			logger.info(String.format("[INTERVAL REDUCE RESULT] at %s ",
					result.getFirst()));
			this.onTraversalDone(result.getSecond());
		}
	}
	 */
	
	
	/**
	 * Return true if your walker wants to reduce each interval separately.
	 * Default is false.
	 * 
	 * If you set this flag, several things will happen.
	 * 
	 * The system will invoke reduceInit() once for each interval being
	 * processed, starting a fresh reduce Reduce will accumulate normally at
	 * each map unit in the interval However, onTraversalDone(reduce) will be
	 * called after each interval is processed. The system will call
	 * onTraversalDone( GenomeLoc -> reduce ), after all reductions are done,
	 * which is overloaded here to call onTraversalDone(reduce) for each
	 * location
	 * 
	 * @return true if your walker wants to reduce each interval separately.
	 */
	@Override
	public boolean isReduceByInterval() {
		return false;
	}



}
