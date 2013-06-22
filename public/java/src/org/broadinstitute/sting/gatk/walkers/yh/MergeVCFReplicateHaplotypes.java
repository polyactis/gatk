/*
 * 2012.5.30
 *  * 	a RODwalker to merge phased calls that are from the same individual. part of the AlignmentToTrioCallPipeline.py workflow.
 * 	an equivalent of phased MergeVCFReplicateGenotypeColumns
 * 
 *  * Note: 1. Samples from VCF that are not in depthFile would be ignored in out.
 * 		2. AC/AF/DP on the INFO column will be updated but QUAL will not.
 * 
 * 
 *  * Examples (--onlyKeepBiAllelicSNP) is optional: java -Xmx2g -jar
 * /home/crocea/script/gatk/dist/GenomeAnalysisTK.jar -T MergeVCFReplicateHaplotypes
 *  -R /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
 *  --variant samtools/Contig0.vcf.gz -o /tmp/contig0_afterMerge.vcf --onlyKeepBiAllelicSNP
 *  --replicateIndividualTag copy
 */
package org.broadinstitute.sting.gatk.walkers.yh;

import java.io.ByteArrayOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;
import java.util.regex.Pattern;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Haplotype;
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
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;


public class MergeVCFReplicateHaplotypes extends RodWalker<Integer, Integer> {
	@ArgumentCollection
	protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

	// 2013.06.18 make this argument optional for the child class CalculateConcordanceAmongReplicates, which does not need it
	protected @Output(doc = "File to which variants should be written", required = false)
	VCFWriter vcfWriter;
	
	@Argument(fullName = "replicateIndividualTag", shortName = "rTag", doc = "tag that separates the true sample ID and the replicate count",
			required = true)
	String replicateIndividualTag = "copy";
	
	// 2011.11.28
	@Argument(fullName = "onlyKeepBiAllelicSNP", doc = "set this to true to remove SNPs with >=3 alleles", required = false)
	Boolean onlyKeepBiAllelicSNP = false;	
	
	// 2012.7.23 debug haplotype distance output file
	private @Argument(fullName = "debugHaplotypeDistanceFname", shortName = "dHDFname", 
			doc = "this file (if given) would contain output of replicate haplotypes to consensus haplotype",
			required = false)
	String debugHaplotypeDistanceFname;
	private FileWriter debugHaplotypeDistanceWriter;
	
	// 2012.7.23 debug output file
	private @Argument(fullName = "debugMajoritySupportFname", shortName = "dMSFname", doc = "this file (if given) would contain output of the support of majority call", required = false)
	String debugMajoritySupportFname;
	private FileWriter debugMajoritySupportWriter;
	
	// key is sample ID, value is a list of minDepth, maxDepth, [minGQ|-1.0].
	Long no_of_total_genotypes = (long) 0;
	Long no_of_included_genotypes = (long) 0;
	Long no_of_retained_loci = (long) 0;
	Long no_of_nonBiAllelicLoci = (long) 0;
	protected TreeSet<String> trueSampleIDSet = new TreeSet<String>();
	public Collection<VariantContext> variantContextCollection = new ArrayList<VariantContext>();	//2012.5.30 to store all genotype data
	public List<DuoConsensusHaplotype> duoConsensusHaplotypeList = new ArrayList<DuoConsensusHaplotype>();	//2012.5.30 each sample has one DuoConsensusHaplotype.

	// 2012.3.30 do keep the original AC, AF, AN
	protected boolean KEEP_ORIGINAL_CHR_COUNTS = true;
	
	Map<String, String> sampleID2trueSampleID = new HashMap<String, String>();
	Map<String, Integer> trueSampleID2Index = new HashMap<String, Integer>();	//2012.5.30 Integer is the index in duoConsensusHaplotypeList
	Map<String, Integer> trueSampleID2Count = new HashMap<String, Integer>();
	Map<String, Diplotype> sampleID2Diplotype = new HashMap<String, Diplotype>();
	int noOfLoci =0;
	
	/*
	 * 2012.7.23 a common function to instruct FileWriter to output a vector of strings
	 */
	public void outputStringVector(FileWriter fw, Vector<String> sv){
		try {
			for (int i = 0; i < sv.size(); i++) {
				fw.append(sv.get(i));
				if (i < sv.size() - 1)
					fw.append("\t");
				else
					fw.append('\n');
			}

		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Error: " + e.getMessage());
		}
	}
	
	/*
	 * 2012.7.25 close the file writer
	 */
	public void closeFileWriter(FileWriter fw){
		if ( fw!= null) {
			try {
				fw.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.err.println("Error: " + e.getMessage());
			}
		}
	}
	
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
			}
			trueSampleID2Count.put(trueSampleID, trueSampleID2Count.get(trueSampleID)+1);
			
			
		}
		
		if (KEEP_ORIGINAL_CHR_COUNTS) {
			headerLines.add(new VCFFormatHeaderLine("AC_Orig", 1,
					VCFHeaderLineType.Integer, "Original AC"));
			headerLines.add(new VCFFormatHeaderLine("AF_Orig", 1,
					VCFHeaderLineType.Float, "Original AF"));
			headerLines.add(new VCFFormatHeaderLine("AN_Orig", 1,
					VCFHeaderLineType.Integer, "Original AN"));
			headerLines.add(new VCFFormatHeaderLine("DP_Orig", 1,
					VCFHeaderLineType.Integer, "Original DP"));
		}
		
		if (vcfWriter!=null){
			//2011-11-11 borrow the header from the input VCF
			vcfWriter.writeHeader(new VCFHeader(headerLines, trueSampleIDSet) );
			//VCFUtils.getHeaderFields(this.getToolkit(),
			//SampleUtils.getUniqueSamplesFromRods(this.getToolkit()) returns empty.
		}
		//2012.7.23
		if (debugMajoritySupportFname != null) {
			try {
				debugMajoritySupportWriter = new FileWriter(debugMajoritySupportFname);
			} catch (IOException e) {
				e.printStackTrace();
				System.err.println("Error: " + e.getMessage());
			}
			//2012.7.23 output the majority vote support & replicate count/family size for each replicate
			Vector<String> majoritySupportDataHeader = new Vector<String>();
			majoritySupportDataHeader.add("#monkeyID");
			majoritySupportDataHeader.add("refStart");
			majoritySupportDataHeader.add("refStop");
			majoritySupportDataHeader.add("noOfLoci");
			majoritySupportDataHeader.add("chromosome");
			majoritySupportDataHeader.add("position");
			majoritySupportDataHeader.add("AF");
			majoritySupportDataHeader.add("hetFraction");
			majoritySupportDataHeader.add("DP");
			majoritySupportDataHeader.add("MQ");
			majoritySupportDataHeader.add("HWE");
			majoritySupportDataHeader.add("VDB");
			majoritySupportDataHeader.add("strandBias");
			majoritySupportDataHeader.add("baseQBias");
			majoritySupportDataHeader.add("mapQBias");
			majoritySupportDataHeader.add("tailDistanceBias");
			
			majoritySupportDataHeader.add("replicateCount");
			majoritySupportDataHeader.add("haplotypeIndex");
			majoritySupportDataHeader.add("base0");
			majoritySupportDataHeader.add("count0");
			majoritySupportDataHeader.add("base1");
			majoritySupportDataHeader.add("count1");
			majoritySupportDataHeader.add("base2");
			majoritySupportDataHeader.add("count2");
			majoritySupportDataHeader.add("base3");
			majoritySupportDataHeader.add("count3");
			outputStringVector(debugMajoritySupportWriter, majoritySupportDataHeader);
				
		
		}
		
		if (debugHaplotypeDistanceFname != null) {
			try {
				debugHaplotypeDistanceWriter = new FileWriter(debugHaplotypeDistanceFname);
			} catch (IOException e) {
				e.printStackTrace();
				System.err.println("Error: " + e.getMessage());
			}
			// output the header for haplotype distance data first
			Vector<String> haplotypeDistanceHeader = new Vector<String>();
			haplotypeDistanceHeader.add("monkeyID_replicateIndex");
			haplotypeDistanceHeader.add("monkeyID");
			haplotypeDistanceHeader.add("replicateCount");
			haplotypeDistanceHeader.add("chromosome");
			haplotypeDistanceHeader.add("start");
			haplotypeDistanceHeader.add("stop");
			haplotypeDistanceHeader.add("refLength");
			haplotypeDistanceHeader.add("noOfLoci");
			haplotypeDistanceHeader.add("replicateSampleIndex");
			haplotypeDistanceHeader.add("haplotype1ToConsensusMinDistance");
			haplotypeDistanceHeader.add("haplotype1ToConsensusMaxDistance");
			haplotypeDistanceHeader.add("haplotype2ToConsensusMinDistance");
			haplotypeDistanceHeader.add("haplotype2ToConsensusMaxDistance");
			outputStringVector(debugHaplotypeDistanceWriter, haplotypeDistanceHeader);
		
		}
	}
	
	public class BaseCount{
		//2012.5.30 class to store the base and the count (frequency)
		private int count;
		private String base;
		public BaseCount(String _base){
			base=_base;
			count = 0;
		}
		
		public void increaseCount(){
			count += 1;
		}
		
		public String getBase(){
			return base;
		}
		public int getCount(){
			return count;
		}
	}
	
	private class BaseCountComparator implements Comparator<BaseCount> {

		public int compare(BaseCount a, BaseCount b) {
			if (a.getCount() < b.getCount())
				return 1;
			if (a.getCount() > b.getCount()) {
				return -1;
			}
			else
				return 0;
		}
	}
	
	public class ConsensusAllele{
		//final PriorityQueue<BaseCount> baseCountQueue = new PriorityQueue<BaseCount>(5, new BaseCountComparator());
		//final Map<String, BaseCount> base2baseCountObject = new HashMap<String, BaseCount>();
		final Map<Byte, Integer> base2Count = new HashMap<Byte, Integer>();
		Byte consensusAllele;
		public Byte returnMostFrequentBase(){
			return consensusAllele;
			//BaseCount baseCountObject = baseCountQueue.peek();
			//return baseCountObject.getBase();
		}
		
		public void addOneBase(String base){
			
			this.addOneBase(Byte.parseByte(base));
		}
		
		public void addOneBase(Byte base){
			if (!base2Count.containsKey(base)){
				base2Count.put(base, 0);
			}
			base2Count.put(base, base2Count.get(base)+1);
			
			int countOfCurrentBase = base2Count.get(base);
			/*
			BaseCount baseCountObject = base2baseCountObject.get(base);
			if (baseCountObject == null){
				baseCountObject = new BaseCount(base);
				base2baseCountObject.put(base, baseCountObject);
			}
			else{
				baseCountQueue.remove(baseCountObject);
			}
			*/
			if (consensusAllele==null){
				consensusAllele = base;
			}
			else{
				int consensusAlleleCount = base2Count.get(consensusAllele);
				if (countOfCurrentBase>consensusAlleleCount){
					consensusAllele = base;
				}
			}
			//baseCountObject.increaseCount();
			//baseCountQueue.add(baseCountObject);
			
		}

	}
	
	
	public class Diplotype{
		Haplotype haplotype1;
		Haplotype haplotype2;
		int noOfLoci;
		ByteArrayOutputStream haplotype1Bases;
		ByteArrayOutputStream haplotype2Bases;
		public byte[] haplotype1ByteArray;	//2013.06.17 more memory efficient than haplotype1, more function than haplotype1Bases
		public byte[] haplotype2ByteArray;
		int currentLocusIndex = 0;
		public Integer refStart;
		public Integer refStop;
		public String chromosome;
		
		public Diplotype(){
			haplotype1Bases = new ByteArrayOutputStream();	//ArrayList<Byte>() uses too much memory
			haplotype2Bases = new ByteArrayOutputStream();	//ArrayList<Byte>();
		}
		public Diplotype(int _noOfLoci){
			noOfLoci = _noOfLoci;
			haplotype1Bases = new ByteArrayOutputStream();	//ArrayList<Byte>();
			haplotype2Bases = new ByteArrayOutputStream();	//ArrayList<Byte>();
		}
		
		public void addGenotypeAtOnePosition(Genotype g, VariantContext vc){
			Allele alleleA = g.getAllele(0);
			Allele alleleB = (g.getAlleles().size() == 2) ? g.getAllele(1) : g.getAllele(0); // hack to deal with no-call genotypes			
			//byte alleleAInByte = alleleA.getBases()[0];
			//byte alleleBInByte = alleleB.getBases()[0];
			haplotype1Bases.write(alleleA.getBaseString().getBytes()[0]);
			haplotype2Bases.write(alleleB.getBaseString().getBytes()[0]);
			
			//#2012.7.25 keep track of the start and stop position of whole haplotype
			if (refStart==null || vc.getStart()<refStart){
				refStart = vc.getStart();
			}
			if (refStop==null || vc.getEnd()>refStop){
				refStop = vc.getEnd();
			}
			if (chromosome==null){
				chromosome = vc.getChr();
			}
			
			currentLocusIndex += 1;
			
		}
		/*
		 * 2013.06.17 return genotype,if sorted=true, then the order of two alleles is sorted:
		 * 	like "AT", instead of 'TA'
		 */
		public byte[] getGenotypeAtOnePosition(int position, boolean sorted){
			
			// check if it's initialized or not
			if (haplotype1ByteArray.length==0){
				initiateHaplotypeByteArray();
			}
			byte[] returnGenotype = new byte[2];
			byte alleleA = haplotype1ByteArray[position-1];
			byte alleleB = haplotype2ByteArray[position-1];
			if (sorted==true){
				if (alleleA>alleleB){
					returnGenotype[0] = alleleB;
					returnGenotype[1] = alleleA;
				}
				else{
					returnGenotype[0] = alleleA;
					returnGenotype[1] = alleleB;
				}
			}
			else{
				returnGenotype[0] = alleleA;
				returnGenotype[1] = alleleB;
			}
			return returnGenotype;
		}
		/*
		 * 2013.06.17 more memory efficient than initiateHaplotypes(), used by CalculateConcordanceAmongReplicates.java
		 * this is needed because haplotype1Bases and haplotype2Bases do not have good methods to access individual allele
		 * 
		 */
		public void initiateHaplotypeByteArray(){
			haplotype1ByteArray = haplotype1Bases.toByteArray();
			haplotype2ByteArray = haplotype2Bases.toByteArray();
			
		}
		/*
		 * this is needed because haplotype1Bases and haplotype2Bases do not have good methods to access individual allele
		 */
		public void initiateHaplotypes(){
			
			//int haplotype1Length = haplotype1Bases.size();
			//byte[] haplotype1InByteArray = ArrayUtils.toPrimitive(haplotype1Bases.toArray(new Byte[0]));
			/*
			byte[] haplotype1InByteArray = haplotype1Bases.toByteArray();
			int haplotype1InByteLength = haplotype1InByteArray.length;
			if (haplotype1Length!=haplotype1InByteLength){
				logger.warn("haplotype1Length " + haplotype1Length + " is not equal to haplotype1InByteLength "+ haplotype1InByteLength);
			}
			*/
			haplotype1 = new Haplotype(haplotype1Bases.toByteArray());
			haplotype2 = new Haplotype(haplotype2Bases.toByteArray());
			//2012.6.12 release some memory
			haplotype1Bases.reset();
			haplotype2Bases.reset();
		}
		
		
	}
	public class DuoConsensusHaplotype{
		public String sampleID;
		public String chromosome;
		public Integer noOfLoci;
		public Integer refStart;
		public Integer refStop;
		ConsensusHaplotypeContainer consensusHaplotype1;
		ConsensusHaplotypeContainer consensusHaplotype2;
		Integer replicateCounter = 0;
		int noOfRandomHaplotypeAssignments = 0;
		int noOfAssignmentsByRelativity = 0;
		Integer noOfReplicates;
		
		public DuoConsensusHaplotype(String _sampleID, int _noOfLoci, String _chromosome){
			sampleID = _sampleID;
			noOfLoci = _noOfLoci;
			chromosome = _chromosome;
			consensusHaplotype1 = new ConsensusHaplotypeContainer(_noOfLoci, 0, _chromosome);
			consensusHaplotype2 = new ConsensusHaplotypeContainer(_noOfLoci, 1, _chromosome);
		}
		

		public DuoConsensusHaplotype(String _sampleID, int _noOfLoci, int _noOfReplicates){
			sampleID = _sampleID;
			noOfLoci = _noOfLoci;
			consensusHaplotype1 = new ConsensusHaplotypeContainer(_noOfLoci, 0);
			consensusHaplotype2 = new ConsensusHaplotypeContainer(_noOfLoci, 1);
			noOfReplicates = _noOfReplicates;
		}
		
		public void addOneDiplotype(Diplotype diplotype){
			Double distance1to1 = consensusHaplotype1.calculateDistanceToConsensusHaplotype(diplotype.haplotype1);
			Double distance2to1 = consensusHaplotype1.calculateDistanceToConsensusHaplotype(diplotype.haplotype2);
			Double distance1to2 = consensusHaplotype2.calculateDistanceToConsensusHaplotype(diplotype.haplotype1);
			Double distance2to2 = consensusHaplotype2.calculateDistanceToConsensusHaplotype(diplotype.haplotype2);
			if (refStart==null){
				refStart = diplotype.refStart;
			}
			if (refStop==null){
				refStop = diplotype.refStop;
			}
			//2012.7.23 output all these distances to debugHaplotypeDistanceWriter if it is not null
			if (debugHaplotypeDistanceWriter !=null){
				if (replicateCounter>0){	//replicateCounter==0 means this diplotype is the 1st entry. meaningless distance
					Vector<String> haplotypeDistanceDataRow = new Vector<String>();
					haplotypeDistanceDataRow.add(sampleID+"_" + replicateCounter.toString());
					haplotypeDistanceDataRow.add(sampleID);
					haplotypeDistanceDataRow.add(noOfReplicates.toString());
					haplotypeDistanceDataRow.add(diplotype.chromosome);
					haplotypeDistanceDataRow.add(diplotype.refStart.toString());
					haplotypeDistanceDataRow.add(diplotype.refStop.toString());
					haplotypeDistanceDataRow.add((new Integer(diplotype.refStop - diplotype.refStart + 1)).toString());
					haplotypeDistanceDataRow.add(noOfLoci.toString());
					haplotypeDistanceDataRow.add(replicateCounter.toString());
					haplotypeDistanceDataRow.add((new Double(Math.min(distance1to1, distance1to2))).toString());
					haplotypeDistanceDataRow.add((new Double(Math.max(distance1to1, distance1to2))).toString());
					haplotypeDistanceDataRow.add((new Double(Math.min(distance2to1, distance2to2))).toString());
					haplotypeDistanceDataRow.add((new Double(Math.max(distance2to1, distance2to2))).toString());
					
					outputStringVector(debugHaplotypeDistanceWriter, haplotypeDistanceDataRow);
					
				}
			}
			double haplotype1DistanceDelta = distance1to2 - distance1to1;	//>0 then haplotype1 is closer to consensus1 than to consensus2
			double haplotype2DistanceDelta = distance2to1 - distance2to2;	//>0 then haplotype2 is closer to consensus2 than to consensus1
			
			replicateCounter += 1;
			if (haplotype1DistanceDelta>0 && haplotype2DistanceDelta>0){
				consensusHaplotype1.addOneHaplotype(diplotype.haplotype1);
				consensusHaplotype2.addOneHaplotype(diplotype.haplotype2);
				
			}
			else if (haplotype1DistanceDelta<0 && haplotype2DistanceDelta<0){
				consensusHaplotype2.addOneHaplotype(diplotype.haplotype1);
				consensusHaplotype1.addOneHaplotype(diplotype.haplotype2);
			}
			else if (haplotype1DistanceDelta<0 && haplotype2DistanceDelta>0){	//both are closer to consensus2
				noOfAssignmentsByRelativity += 1;
				if (-haplotype1DistanceDelta<haplotype2DistanceDelta){	//haplotype1 is closer to consensus1 than haplotype1.
					consensusHaplotype1.addOneHaplotype(diplotype.haplotype1);
					consensusHaplotype2.addOneHaplotype(diplotype.haplotype2);
				}
				else{
					consensusHaplotype2.addOneHaplotype(diplotype.haplotype1);
					consensusHaplotype1.addOneHaplotype(diplotype.haplotype2);
				}
			}
			else if (haplotype1DistanceDelta>0 && haplotype2DistanceDelta<0){	//both are closer to consensus1
				noOfAssignmentsByRelativity += 1;
				if (haplotype1DistanceDelta<-haplotype2DistanceDelta){	//haplotype2 is closer to consensus1 than haplotype1.
					consensusHaplotype2.addOneHaplotype(diplotype.haplotype1);
					consensusHaplotype1.addOneHaplotype(diplotype.haplotype2);
				}
				else{
					consensusHaplotype1.addOneHaplotype(diplotype.haplotype1);
					consensusHaplotype2.addOneHaplotype(diplotype.haplotype2);
				}
			}
			else if ((haplotype1DistanceDelta>0 && haplotype2DistanceDelta==0) || (haplotype1DistanceDelta==0 && haplotype2DistanceDelta>0)){
				consensusHaplotype1.addOneHaplotype(diplotype.haplotype1);
				consensusHaplotype2.addOneHaplotype(diplotype.haplotype2);
			}
			else if ((haplotype1DistanceDelta<0 && haplotype2DistanceDelta==0) || (haplotype1DistanceDelta==0 && haplotype2DistanceDelta<0)){
				consensusHaplotype2.addOneHaplotype(diplotype.haplotype1);
				consensusHaplotype1.addOneHaplotype(diplotype.haplotype2);
			}
			else{	//both are 0
				consensusHaplotype1.addOneHaplotype(diplotype.haplotype1);
				consensusHaplotype2.addOneHaplotype(diplotype.haplotype2);
				noOfRandomHaplotypeAssignments += 1;
			}
		}
		
		public Genotype getGenotypeAtOnePosition(int i, String refString, Genotype oldGenotype){
			
			String alleleA = new String(new byte[] {consensusHaplotype1.getOneBaseFromConsensusHaplotype(i)});
			String alleleB = new String(new byte[] {consensusHaplotype2.getOneBaseFromConsensusHaplotype(i)});
			
			//Allele originalAlleleB = (g.getAlleles().size() == 2) ? g.getAllele(1) : g.getAllele(0); // hack to deal with no-call genotypes
			boolean genotypeIsPhased = true;
			List<Allele> alleles = new ArrayList<Allele>();
			Allele bglAlleleA, bglAlleleB;

			if (alleleA.matches(refString))
				bglAlleleA = Allele.create(alleleA, true);
			else
				bglAlleleA = Allele.create(alleleA, false);

			if (alleleB.matches(refString))
				bglAlleleB = Allele.create(alleleB, true);
			else
				bglAlleleB = Allele.create(alleleB, false);
			alleles.add(bglAlleleA);
			alleles.add(bglAlleleB);
			
			HashMap<String, Object> originalAttributes = new HashMap<String, Object>(oldGenotype.getAttributes());
			Set<String> filters = new LinkedHashSet<String>(oldGenotype.getFilters());
			double genotypeQuality = oldGenotype.getLog10PError();
			
			Genotype genotype = new Genotype(sampleID, alleles, genotypeQuality, filters, originalAttributes , genotypeIsPhased);
			return genotype;
			
		}
		/*
		 * 2012.7.23 
		 */
		public ConsensusAllele getConsensusAlleleAtOnePosition(int i, int haplotypeIndex){
			if (haplotypeIndex==0){
				return consensusHaplotype1.consensusAlleleList.get(i);
			}
			else{
				return consensusHaplotype2.consensusAlleleList.get(i);
			}
		}
		public void reportConsensusStat(){
			logger.info(replicateCounter + " replicates for sample " + sampleID);
			logger.info("Haplotypes of "+ noOfRandomHaplotypeAssignments + " replicates were randomly assigned.");
			logger.info("Haplotypes of "+ noOfAssignmentsByRelativity + " replicates were assigned by relativity.");
		}
		
		
	}
	/*
	 * 2012.5.31 to be unit of DuoConsensusHaplotype
	 */
	public class ConsensusHaplotypeContainer{
		public String chromosome;
		public Haplotype consensusHaplotype;
		public List<ConsensusAllele> consensusAlleleList;	//one for each locus
		public int noOfLoci;
		public Integer phaseIndex;	//2012.7.24 which chromosome it is on.
		public Integer noOfReplicates=0;	//2012.7.25
		public ConsensusHaplotypeContainer(int _noOfLoci, int _phaseIndex, String _chromosome){
			noOfLoci = _noOfLoci;
			phaseIndex = _phaseIndex;
			chromosome = _chromosome;
		}

		public ConsensusHaplotypeContainer(int _noOfLoci, int _phaseIndex){
			noOfLoci = _noOfLoci;
			phaseIndex = _phaseIndex;
			
		}
		
		/*
		 * 2012.7.25
		 */
		public void initConsensusAlleleList(Haplotype haplotype){
			consensusAlleleList = new ArrayList<ConsensusAllele>();	//one for each locus
			int _noOfLoci = haplotype.getBases().length;
			for (int i=0; i<_noOfLoci; i++){
				ConsensusAllele consensusAllele = new ConsensusAllele();
				consensusAllele.addOneBase(haplotype.getBases()[i]);
				consensusAlleleList.add(consensusAllele);
			}
			
		}
		public double calculateDistanceToConsensusHaplotype(Haplotype haplotype){
			double noOfMismatches = 0.0;
			if (consensusHaplotype !=null){
				for (int i=0; i<haplotype.getBases().length;i++){
					if (haplotype.getBases()[i]!=consensusHaplotype.getBases()[i]){
						noOfMismatches += 1;
					}
				}
				return noOfMismatches/(double)noOfLoci;
			}
			else{
				return 0.0;
			}
			
		}
		
		private void updateConsensusHaplotype(){
			if (consensusAlleleList!=null){
				final byte[] haplotypeBases = new byte[noOfLoci];
				for (int i=0; i<consensusAlleleList.size(); i++){
					haplotypeBases[i] = consensusAlleleList.get(i).returnMostFrequentBase();
				}
				
				consensusHaplotype = new Haplotype(haplotypeBases);
			}
		}
		/*
		 * 2012.5.30
		 */		
		public void addOneHaplotype(Haplotype haplotype){
			int _noOfLoci = haplotype.getBases().length;
			if (noOfReplicates>0){	//2012.7.25 second replicate comes in and use consensusAlleleList to hold data
				initConsensusAlleleList(consensusHaplotype);	//initialize the consensusAlleleList with consensusHaplotype (1st haplotype) 
				
				for (int i=0; i<_noOfLoci; i++){
					ConsensusAllele consensusAllele = consensusAlleleList.get(i);
					consensusAllele.addOneBase(haplotype.getBases()[i]);
				}
				updateConsensusHaplotype();
			}
			else{	//2012.7.25 1st haplotype. directly push its bases into consensusHaplotype.
				final byte[] haplotypeBases = new byte[noOfLoci];
				for (int i=0; i<_noOfLoci; i++){
					haplotypeBases[i] = haplotype.getBases()[i];
				}
				consensusHaplotype = new Haplotype(haplotypeBases);
			}
			noOfReplicates ++;	//2012.7.25
			
		}
		
		public Byte getOneBaseFromConsensusHaplotype(int i){
			return consensusHaplotype.getBases()[i];
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
			noOfLoci += 1;
			variantContextCollection.add(vc);
			GenotypesContext thisGC = vc.getGenotypes();
			for (Genotype genotype: thisGC) {
				no_of_total_genotypes ++;
				String sampleID = genotype.getSampleName();
				if (!sampleID2Diplotype.containsKey(sampleID)){
					sampleID2Diplotype.put(sampleID, new Diplotype());
				}
				Diplotype diplotype = sampleID2Diplotype.get(sampleID);
				diplotype.addGenotypeAtOnePosition(genotype, vc);
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

	
	public void onTraversalDone(Integer result) {
		int i=0;
		
		if ( duoConsensusHaplotypeList.isEmpty()){
			Iterator<String> itr = trueSampleIDSet.iterator();
			i=0;
			while (itr.hasNext()){
				String sampleID = (String) itr.next();
				duoConsensusHaplotypeList.add(new DuoConsensusHaplotype(sampleID, noOfLoci, trueSampleID2Count.get(sampleID)));
				trueSampleID2Index.put(sampleID, i);
				i += 1;
			}
		}
		// 2012.5.30 initiate each haplotype and add them to duoConsensusHaplotype
		Set<String> sampleIDSet = sampleID2Diplotype.keySet();
		for (String sampleID: sampleIDSet){
			sampleID2Diplotype.get(sampleID).initiateHaplotypes();
			String trueSampleID = sampleID2trueSampleID.get(sampleID);
			int dataStructureIndex = trueSampleID2Index.get(trueSampleID);
			DuoConsensusHaplotype duoConsensusHaplotype = duoConsensusHaplotypeList.get(dataStructureIndex);
			duoConsensusHaplotype.addOneDiplotype(sampleID2Diplotype.get(sampleID));
			//sampleID2Diplotype.remove(sampleID);	//2012.7.25 it'll result in java.util.ConcurrentModificationException
		}
		//2012.6.12 release memory
		sampleID2Diplotype.clear();
		
		
		//2012.5.30 finally to output
		i = -1;	//2012.5.31 starting from -1 because "+=" is executed right after for loop to make sure it's incremented even when "continue" is executed
		for (VariantContext vc: variantContextCollection){
			// loop through each locus (row)
			i += 1;
			if (onlyKeepBiAllelicSNP==true){
				//2012.5.31
				if (vc.getAlternateAlleles().size()!=1){
					no_of_nonBiAllelicLoci ++;
					continue;
				}
			}
			
			GenotypesContext thisGC = vc.getGenotypes();
			//2012.7.23 trueSampleID2OneGenotype will be template for the new genotype.
			Map<String, Genotype> trueSampleID2OneGenotype = new HashMap<String, Genotype>();
			for (Genotype genotype: thisGC) {
				no_of_total_genotypes ++;
				String sampleID = genotype.getSampleName();
				String trueSampleID = sampleID2trueSampleID.get(sampleID);
				if (!trueSampleID2OneGenotype.containsKey(trueSampleID)){
					trueSampleID2OneGenotype.put(trueSampleID, genotype);
				}
			}
			
			GenotypesContext newGC = GenotypesContext.create();
			// newGC.clear();
			String refString = vc.getReference().getDisplayString();
			if (refString.length() == 0) // ref was null
				refString = Allele.NULL_ALLELE_STRING;
			int depth = 0;
			for (Map.Entry<String, Integer> item: trueSampleID2Index.entrySet()){
				// loop through each sample (column or individual)
				String trueSampleID = item.getKey();
				int duoConsensusHaplotypeIndex = item.getValue();
				DuoConsensusHaplotype duoConsensusHaplotype = duoConsensusHaplotypeList.get(duoConsensusHaplotypeIndex);
				//if (i==0){	//#first locus, report 
				//	duoConsensusHaplotype.reportConsensusStat();
				//}
				Genotype oldGenotype = trueSampleID2OneGenotype.get(trueSampleID);
				Genotype genotype = duoConsensusHaplotype.getGenotypeAtOnePosition(i, refString, oldGenotype);
				newGC.add(genotype);
				Integer dp = getGenotypeDepth(genotype);
				depth += dp;
				if (debugMajoritySupportWriter!=null && duoConsensusHaplotype.replicateCounter>1){
					List<ConsensusAllele> consensusAlleleList = new ArrayList<ConsensusAllele>();	//one for each locus
					
					consensusAlleleList.add(duoConsensusHaplotype.getConsensusAlleleAtOnePosition(i, 0));
					consensusAlleleList.add(duoConsensusHaplotype.getConsensusAlleleAtOnePosition(i, 1));
					for (Integer m=0;m<consensusAlleleList.size(); m++){
						ConsensusAllele haplotypeConsensusAllele = consensusAlleleList.get(m);
						if (haplotypeConsensusAllele.base2Count.size()>1){
							Vector<String> majoritySupportDataRow  = new Vector<String>();
							majoritySupportDataRow.add(trueSampleID);
							majoritySupportDataRow.add(duoConsensusHaplotype.refStart.toString());
							majoritySupportDataRow.add(duoConsensusHaplotype.refStop.toString());
							majoritySupportDataRow.add(duoConsensusHaplotype.noOfLoci.toString());	//loci-range
							majoritySupportDataRow.add(vc.getChr());
							majoritySupportDataRow.add(""+vc.getStart());	//adding an empty string converts the int to string
							
							majoritySupportDataRow.add(getVariantAF(vc).toString());

							Double hetFraction = vc.getHetCount()
									/ (double) vc.getNSamples();
							majoritySupportDataRow.add(hetFraction.toString());

							Integer DP = vc.getAttributeAsInt("DP", -1);
							majoritySupportDataRow.add(DP.toString());

							Integer mq = vc.getAttributeAsInt("MQ", -1);
							majoritySupportDataRow.add(mq.toString());

							majoritySupportDataRow.add(vc.getAttributeAsString("HWE", "-1"));

							majoritySupportDataRow.add(vc.getAttributeAsString("VDB", "-1"));

							String pv4 = vc.getAttributeAsString("PV4", "");
							String[] pv4List = new String[] { "-1", "-1", "-1", "-1" };
							if (pv4.length() > 0) {
								pv4List = pv4.substring(1, pv4.length() - 1).split(",");
							}
							for (int n = 0; n < pv4List.length; n++) {
								majoritySupportDataRow.add(pv4List[n]);
							}
							
							majoritySupportDataRow.add(duoConsensusHaplotype.replicateCounter.toString());
							majoritySupportDataRow.add(m.toString());	//haplotype index
							for (Map.Entry<Byte, Integer> base_count_item: haplotypeConsensusAllele.base2Count.entrySet()){
								majoritySupportDataRow.add(base_count_item.getKey().toString());
								Double fraction = base_count_item.getValue().doubleValue()/duoConsensusHaplotype.replicateCounter.doubleValue();
								majoritySupportDataRow.add(fraction.toString());
							}
							outputStringVector(debugMajoritySupportWriter, majoritySupportDataRow);
						}
					}
					// do the same thing for the haplotype1ConsensusAllele
				}
			}
			// 2012.6.1 vc.getAlleles() includes only one alternative alleles (tri-alleles). cuz failure when a 3rd allele shows up in newGC.
			List<Allele> allAlleles = new ArrayList<Allele>();
			allAlleles.add(vc.getReference());
			allAlleles.addAll(vc.getAlternateAlleles());
			
			VariantContext sub = vc.subContextFromSamples(trueSampleIDSet, allAlleles);
			VariantContextBuilder builder = new VariantContextBuilder(sub);
			builder.genotypes(newGC);
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
				if (sub.hasAttribute(VCFConstants.DEPTH_KEY))
					builder.attribute("DP_Orig",
							sub.getAttribute(VCFConstants.DEPTH_KEY));
			}
			builder.attribute("DP", depth);
			VariantContextUtils.calculateChromosomeCounts(builder, false);
			VariantContext newVC = builder.make();
			
			
			
			if (newVC.getGenotypes().size()>0){	//output only when not all are missing
				no_of_retained_loci ++;
				if (vcfWriter!=null){
					vcfWriter.add(newVC);
				}
			}
		}
		
		//2012.7.24 close the file, otherwise it might end up with some missing end (buffer not flushed to disk)
		closeFileWriter(debugHaplotypeDistanceWriter);
		closeFileWriter(debugMajoritySupportWriter);
		
		// 2012.7.24 don't try to close vcfWriter (it'll close itself somewhere after)
		logger.info(result + " loci processed.");
		logger.info(noOfLoci + " loci in each haplotype.");
		logger.info(no_of_retained_loci + " loci retained.");
		logger.info("Number of total genotypes: " + no_of_total_genotypes);
		logger.info("Number of outputted genotypes: "
				+ no_of_included_genotypes);
		logger.info("Number of non-bi-allelic loci after GQ and Depth filter: "
				+ no_of_nonBiAllelicLoci);
	}
	
}
