package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

import aligner.SingleStateAlignerFlat2;
import dna.AminoAcid;
import fileIO.FileFormat;
import gff.GffLine;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;
import structures.IntList;

/**
 * This class is designed to store kmer frequencies related to gene
 * starts, stops, and interiors.  It can be loaded from a pgm file.
 * 
 * It's possible to use multiple GeneModels; for example, one for
 * each of several GC ranges or clades.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 *
 */
public class GeneModel extends ProkObject {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public GeneModel(boolean fill){
		if(fill){
			fillContainers();
		}
	}
	
	void fillContainers(){
		statsCDS.setInner(kInnerCDS, 3);
		statsCDS.setStart(kStartCDS, startFrames, startLeftOffset);
		statsCDS.setStop(kStopCDS, stopFrames, stopLeftOffset);

		for(int i=0; i<rnaContainers.length; i++){
			StatsContainer sc=rnaContainers[i];
			sc.setInner(kInnerRNA, 1);
		}

		statstRNA.setStart(kStartRNA, 14, 4);
		statstRNA.setStop(kStopRNA, 14, 6);

		stats16S.setStart(kStartRNA, 20, 7);
		stats16S.setStop(kStopRNA, 12, 16);

		stats23S.setStart(kStartRNA, 17, 3);
		stats23S.setStop(kStopRNA, 15, 12);

		stats5S.setStart(kStartRNA, 20, 5);
		stats5S.setStop(kStopRNA, 15, 5);

		stats18S.setStart(kStartRNA, 20, 7);//TODO: 18S bounds are untested and should be empirically determined
		stats18S.setStop(kStopRNA, 12, 16);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean process(String genomeFname, String gffFname){
		FileFormat fnaFile=FileFormat.testInput(genomeFname, FileFormat.FA, null, true, true);
		FileFormat gffFile=FileFormat.testInput(gffFname, FileFormat.GFF, null, true, true);
		
		if(fnaFile==null || gffFile==null){
			errorState=true;
			return true;
		}
		
		ArrayList<Read> reads=ReadInputStream.toReads(fnaFile, maxReads);
		ArrayList<GffLine>[] allGffLines=GffLine.loadGffFileByType(gffFile, "CDS,rRNA,tRNA", true);
		ArrayList<GffLine> cds=allGffLines[0];
		ArrayList<GffLine> rrna=allGffLines[1];
		ArrayList<GffLine> trna=allGffLines[2];

		return process(reads, cds, rrna, trna);
	}
	
	public boolean process(ArrayList<Read> genome, 
			ArrayList<GffLine> cds, ArrayList<GffLine> rrna, ArrayList<GffLine> trna){
		numFiles++;
		
		if(genome==null || cds==null){
			errorState=true;
			return true;
		}
		filesProcessed++;
		
		ArrayList<ScafData> scafList;
		{//Scoped to save memory
			readsProcessed+=genome.size();
			scafList=new ArrayList<ScafData>(genome.size());
			for(Read r : genome){
				basesProcessed+=r.length();
				scafList.add(new ScafData(r));
			}
		}
		{//Scoped to save memory
			genesProcessed+=cds.size();
			genesProcessed+=(rrna==null ? 0 : rrna.size());
			genesProcessed+=(trna==null ? 0 : trna.size());
			
			HashMap<String, ScafData> scafMap=makeScafMap(scafList);
			fillScafDataCDS(cds, scafMap);
			fillScafDataRNA(rrna, scafMap);
			fillScafDataRNA(trna, scafMap);
		}
		
		countBases(scafList);
		if(PROCESS_PLUS_STRAND){
			processStrand(scafList, Shared.PLUS);
		}
		if(PROCESS_MINUS_STRAND){
			for(ScafData sd : scafList){
				sd.clear();
				sd.reverseComplement();
			}
			processStrand(scafList, Shared.MINUS);
			for(ScafData sd : scafList){
				sd.clear();
				sd.reverseComplement();
			}
		}
		return false;
	}
	
	public static boolean ADD_CDS_ONLY=false;//Sometimes merging in a custom model messes up non-CDS numbers
	public void add(GeneModel pgm){
		if(ADD_CDS_ONLY && readsProcessed>0) {
//			System.err.println("Adding CDS; readsProcessed="+readsProcessed+", ADD_CDS_ONLY="+ADD_CDS_ONLY);
			statsCDS.add(pgm.statsCDS);
		}else {
//			System.err.println("Adding all; readsProcessed="+readsProcessed+", ADD_CDS_ONLY="+ADD_CDS_ONLY);
			for(int i=0; i<allContainers.length; i++){
//				System.err.println("merging "+allContainers[i].name);
				allContainers[i].add(pgm.allContainers[i]);
			}
		}
		
		readsProcessed+=pgm.readsProcessed;
		basesProcessed+=pgm.basesProcessed;
		genesProcessed+=pgm.genesProcessed;
		filesProcessed+=pgm.filesProcessed;

//		geneStartsProcessed+=pgm.geneStartsProcessed;
//		tRNAProcessed+=pgm.tRNAProcessed;
//		r16SProcessed+=pgm.r16SProcessed;
//		r23SProcessed+=pgm.r23SProcessed;
//		r5SProcessed+=pgm.r5SProcessed;
//		r18SProcessed+=pgm.r18SProcessed;
		
//		fnames.addAll(pgm.fnames);
		numFiles+=pgm.numFiles;
		taxIds.addAll(pgm.taxIds);
		Tools.add(baseCounts, pgm.baseCounts);
	}
	
	public void multiplyBy(double mult) {
		for(int i=0; i<allContainers.length; i++){
			allContainers[i].multiplyBy(mult);
		}
		
		readsProcessed=Math.round(readsProcessed*mult);
		basesProcessed=Math.round(basesProcessed*mult);
		genesProcessed=Math.round(genesProcessed*mult);
		filesProcessed=Math.round(filesProcessed*mult);
		
		for(int i=0; i<baseCounts.length; i++){
			baseCounts[i]=Math.round(baseCounts[i]*mult);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public float gc(){
		long a=baseCounts[0];
		long c=baseCounts[1];
		long g=baseCounts[2];
		long t=baseCounts[3];
		return (float)((g+c)/Tools.max(1.0, a+t+g+c));
	}
	
	HashMap<String, ScafData> makeScafMap(ArrayList<ScafData> scafList){
		HashMap<String, ScafData> scafMap=new HashMap<String, ScafData>(scafList.size()*3);
		for(ScafData sd : scafList){scafMap.put(sd.name, sd);}
		for(ScafData sd : scafList){
			String name=sd.name;
			int idx=name.indexOf(' ');
			if(idx>=0){
				String prefix=name.substring(0, idx);
				if(scafMap.containsKey(prefix)){
					assert(false) : "Duplicate degenerate name: '"+name+"', '"+prefix+"'";
				}else{
					scafMap.put(prefix, sd);
				}
			}
		}
		return scafMap;
	}
	
	public void fillScafDataCDS(ArrayList<GffLine> cdsLines, HashMap<String, ScafData> scafMap){
		if(!callCDS){return;}
		for(GffLine gline : cdsLines){
			ScafData sd=scafMap.get(gline.seqid);
			assert(sd!=null) : "Can't find scaffold for GffLine "+gline.seqid;
			sd.addCDS(gline);
		}
	}
	
	public void fillScafDataRNA(ArrayList<GffLine> rnaLines, HashMap<String, ScafData> scafMap){
		for(GffLine gline : rnaLines){
			ScafData sd=scafMap.get(gline.seqid);
			assert(sd!=null) : "Can't find scaffold for GffLine "+gline.seqid;
			if(processType(gline.prokType())){
				sd.addRNA(gline);
			}
		}
	}
	
	public void processStrand(ArrayList<ScafData> scafList, int strand){
		for(ScafData sd : scafList){
			processCDS(sd, strand);
			processRNA(sd, strand);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	private void countBases(ArrayList<ScafData> scafList){
		for(ScafData sd : scafList){
			countBases(sd.bases);
		}
	}
	
	private void countBases(byte[] bases){
		for(byte b : bases){
			int x=AminoAcid.baseToNumberACGTother[b];
			baseCounts[x]++;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Finding Codons        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static void findStopCodons(byte[] bases, IntList list, BitSet valid){
		final int k=3;
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					int point=i;//End of the stop codon
					if(isStopCodon(kmer) && !valid.get(point)){
						list.add(point);
						valid.set(point);
					}
				}
			}else{len=0;}
		}
		
		for(int i=50; i<bases.length-3; i+=2000){//Add some non-canonical sites, aka noise
			if(!valid.get(i)){
				list.add(i);
			}
		}
	}
	
	private static void findStartCodons(byte[] bases, IntList list, BitSet valid){
		final int k=3;
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					int point=i-k+1;//Start of the start codon
					if(isStartCodon(kmer) && !valid.get(point)){
						list.add(point);
						valid.set(point);
					}
				}
			}else{len=0;}
		}
		
		for(int i=50; i<bases.length-3; i+=2000){//Add some non-canonical sites, aka noise
			if(!valid.get(i)){
				list.add(i);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Processing GffLines      ----------------*/
	/*--------------------------------------------------------------*/
	
	private static void processGene(GffLine gline, ScafData sd){
		if(gline.length()<2){return;}
		final int strand=gline.strand;
		assert(strand==sd.strand());
		final byte[] frames=sd.frames;
		int start=gline.start-1, stop=gline.stop-1;
		if(start<0 || stop>=sd.length()){return;}
//		assert(start<stop) : gline; //Not always true for euks...
		if(strand==Shared.MINUS){
			int x=sd.length()-start-1;
			int y=sd.length()-stop-1;
			start=y;
			stop=x;

//			String a=new String(sd.bases, start, 3);
//			String b=new String(sd.bases, stop-2, 3);
////			assert(false) : start+", "+stop+"\n"+gline+"\n"+new String(sd.bases, start, 3)+", "+new String(sd.bases, stop-2, 3);
//			outstream.println(a+", "+b+", "+start+", "+stop);
		}
		assert(start>=0) : gline.toString()+"\n"+sd.length()+"\n"+sd.name;
		markFrames(start, stop, frames, kInnerCDS);
		sd.starts.add(start);
		sd.stops.add(stop);
//		assert(gline.start!=337) : gline+"\n"+start+", "+stop;
	}
	
	private boolean processRnaLine(final GffLine gline, final ScafData sd, final int type){
		final int strand=gline.strand;
		assert(strand==sd.strand());
		final byte[] frames=sd.frames;
		int start=gline.start-1, stop=gline.stop-1;
		if(start<0 || stop>=sd.length()){return false;}
		assert(start<=stop);
		if(strand==Shared.MINUS){
			int x=sd.length()-start-1;
			int y=sd.length()-stop-1;
			start=y;
			stop=x;
		}
		
		if(AnalyzeGenes.alignRibo){
//			byte[] seq=sd.fetch(start, stop);
			Read[] consensusReads=ProkObject.consensusReads(type);
			byte[] universal=(consensusReads!=null && consensusReads.length>0 ? consensusReads[0].bases : null);
			float minIdentity=ProkObject.minID(type);
			if(universal!=null){
				int[] coords=KillSwitch.allocInt1D(2);
				final int a=Tools.max(0, start-(AnalyzeGenes.adjustEndpoints ? 200 : 50));
				final int b=Tools.min(sd.bases.length-1, stop+(AnalyzeGenes.adjustEndpoints ? 200 : 50));
				float id1=align(universal, sd.bases, a, b, minIdentity, coords);
				final int rstart=coords[0], rstop=coords[1];
//				assert(false) : rstart+", "+rstop+", "+(rstop-rstart+1)+", "+start+", "+stop;
				if(id1<minIdentity){
//					System.err.println("Low identity: "+String.format("%.2s", 100*id1));
					return false;
				}else{
//					System.err.println("Good identity: "+String.format("%.2s", 100*id1));
				}
				if(AnalyzeGenes.adjustEndpoints){
					int startSlop=startSlop(type);
					int stopSlop=stopSlop(type);
					if(Tools.absdif(start, rstart)>startSlop){
//						System.err.println("rstart:\t"+start+" -> "+rstart);
						start=rstart;
					}
					if(Tools.absdif(stop, rstop)>stopSlop){
//						System.err.println("rstop: \t"+stop+" -> "+rstop);
						stop=rstop;
					}
				}
			}
		}
		
		StatsContainer sc=allContainers[type];
		sc.start.processPoint(sd.bases, start, 1);
		sc.stop.processPoint(sd.bases, stop, 1);
		assert(sc!=statsCDS);
		
		assert(start>=0) : gline.toString()+"\n"+sd.length()+"\n"+sd.name;
		final byte flag=typeToFlag(type);
		for(int i=start+kInnerRNA-1; i<=stop; i++){
			frames[i]|=flag;
		}
		return true;
	}
	
	private float align(byte[] query, byte[] ref, int start, int stop, float minIdentity, int[] coords){
//		final int a=0, b=ref.length-1;
		SingleStateAlignerFlat2 ssa=GeneCaller.getSSA();
		final int minScore=ssa.minScoreByIdentity(query.length, minIdentity);
		int[] max=ssa.fillUnlimited(query, ref, start, stop, minScore);
		if(max==null){return 0;}
		
		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
//		final int maxScore=max[3];
//		final int maxStart=max[4];
		
		//returns {score, bestRefStart, bestRefStop} 
		//padded: {score, bestRefStart, bestRefStop, padLeft, padRight};
		final int[] score=ssa.score(query, ref, start, stop, rows, maxCol, maxState);
		final int rstart=Tools.max(score[1], 0);
		final int rstop=Tools.min(score[2], ref.length-1);
		if(coords!=null){
			coords[0]=rstart;
			coords[1]=rstop;
		}
		final float id=ssa.tracebackIdentity(query, ref, start, stop, rows, maxCol, maxState, null);
		return id;
	}
	
	/** 
	 * Each frame byte has a bit marked for valid coding frames.
	 * For example, if frames[23]=0b100, then base 23 is the last base in a kmer starting at the 3rd base in a codon.
	 * If frames[23]=0, then no coding kmer end at that location on this strand.
	 * @param start
	 * @param stop
	 * @param frames
	 * @param k
	 */
	private static void markFrames(int start, int stop, byte[] frames, int k){
		assert(start<=stop) : start+", "+stop;
		for(int i=start+k-1, frameBit=(1<<((k-1)%3)), max=Tools.min(stop-3, frames.length-1); i<=max; i++){
			frames[i]=(byte)(frames[i]|frameBit);
			frameBit<<=1;
			if(frameBit>4){frameBit=1;}
		}
//		assert(false) : Arrays.toString(Arrays.copyOfRange(frames, start, start+20))+"\n"+start; //This is correct
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Counting Kmers        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void processCDS(ScafData sd, int strand){
		if(!callCDS){return;}
		ArrayList<GffLine> glines=sd.cdsLines[strand];
		for(GffLine gline : glines){
			assert(gline.strand==strand);
			processGene(gline, sd);
			statsCDS.addLength(gline.length());
		}
		
		statsCDS.inner.processCDSFrames(sd.bases, sd.frames);
		BitSet startSet=processEnds(sd.bases, statsCDS.start, sd.starts, 1);
		BitSet stopSet=processEnds(sd.bases, statsCDS.stop, sd.stops, 1);
//		outstream.println("Processed "+sd.starts.size+" valid starts and "+sd.stops.size+" stops.");
		sd.clear();
		findStartCodons(sd.bases, sd.starts, startSet);
		findStopCodons(sd.bases, sd.stops, stopSet);
//		outstream.println("Found "+sd.starts.size+" invalid starts and "+sd.stops.size+" stops.");
		processEnds(sd.bases, statsCDS.start, sd.starts, 0);
		processEnds(sd.bases, statsCDS.stop, sd.stops, 0);
	}
	
	private static int getGlineType(GffLine gline, ScafData sd){
		if(!gline.inbounds(sd.length()) || gline.partial()){return -1;}

		final int length=gline.length();
		final int type=gline.prokType();
		if(type<0){
			return type;
		}else if(type==CDS){
			return type;
		}else if(type==tRNA && length>=40 && length<=120){
			return type;
		}else if(type==r16S && length>=1440 && length<=1620){
			return type;
		}else if(type==r23S && length>=2720 && length<=3170){
			return type;
		}else if(type==r5S && length>=90 && length<=150){
			return type;
		}else if(type==r18S && length>=1400 && length<=2000){ //TODO: Check length range
			return type;
		}
		return -1;
	}
	
	private void processRNA(ScafData sd, int strand){
		sd.clear();
		ArrayList<GffLine> lines=sd.rnaLines[strand];
		for(GffLine gline : lines){
			assert(gline.strand==strand);
			final int type=getGlineType(gline, sd);
			if(type>0){
				StatsContainer sc=allContainers[type];
				sc.addLength(gline.length());
				processRnaLine(gline, sd, type);
			}
		}
		processRnaInner(sd);
		processRnaEnds(sd);
	}
	
	void processRnaInner(ScafData sd){
		byte[] bases=sd.bases;
		byte[] frames=sd.frames;
		final int k=kInnerRNA;//TODO: Note! This is linked to a single static variable for all RNAs.
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					int vf=frames[i];
					for(int type=0; type<5; type++){
						int valid=vf&1;
						rnaContainers[type].inner.add(kmer, 0, valid);
						vf=(vf>>1);
					}
				}
			}else{len=0;}
		}
	}
	
	void processRnaEnds(ScafData sd){
		byte[] bases=sd.bases;

		final int k=stats16S.start.k;
		final int kMax=stats16S.start.kMax;
		final int mask=stats16S.start.mask;
		final long[] counts=new long[kMax];//TODO: Slow
		
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					counts[kmer]++;
				}
			}else{len=0;}
		}
		for(StatsContainer sc : rnaContainers){
			FrameStats fs=sc.start;
			for(long[] array : fs.countsFalse){
				Tools.add(array, counts);
			}
			fs=sc.stop;
			for(long[] array : fs.countsFalse){
				Tools.add(array, counts);
			}
		}
	}
	
	private static BitSet processEnds(byte[] bases, FrameStats stats, IntList list, int valid){
		BitSet points=new BitSet(bases.length);
		for(int i=0; i<list.size; i++){
			int point=list.get(i);
			stats.processPoint(bases, list.get(i), valid);
			points.set(point);
		}
		return points;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Scoring            ----------------*/
	/*--------------------------------------------------------------*/
	
//	//Assumes bases are in the correct strand
//	public float calcStartScore(int start, int stop, byte[] bases){
//		float f=scorePoint(start, bases, startStats);
////		float ss=scoreStart2(start, bases, stop, innerKmerStats);
////		if(ss>0){f=(f+0.0005f*ss);} //Does not seem to help; needs more study.
//		return f;
//	}
//	
//	//Assumes bases are in the correct strand
//	public float calcStopScore(int stop, byte[] bases){
//		float f=scorePoint(stop, bases, stopStats);
//		return f;
//	}
//	
//	//Assumes bases are in the correct strand
//	public float calcRnaStartScore(int start, int stop, byte[] bases){
//		float f=scorePoint(start, bases, rrnaStartStats);
//		return f;
//	}
//	
//	//Assumes bases are in the correct strand
//	public float calcRnaStopScore(int stop, byte[] bases){
//		float f=scorePoint(stop, bases, rrnaStopStats);
//		return f;
//	}
	
//	public static float calcKmerScore(int start, int stop, int startFrame, byte[] bases, FrameStats stats){
//
//		assert(stats.frames==3);
//		final int k=stats.k;
//		final int mask=~((-1)<<(2*k));
//		
//		int kmer=0;
//		int len=0;
//		float score=0;
//		int numKmers=0;
//		
//		for(int pos=start, currentFrame=startFrame; pos<stop; pos++){
//			final byte b=bases[pos];
//			final int x=AminoAcid.baseToNumber[b];
//
//			if(x>=0){
//				kmer=((kmer<<2)|x)&mask;
//				len++;
//				if(len>=k){
//					float prob=stats.probs[currentFrame][kmer];
//					float dif=prob-0.99f;//Prob above 1 is more likely than average
//					score+=dif;
//					numKmers++;
//				}
//			}else{
//				len=0;
//				kmer=0;
//			}
//
//			currentFrame++;
//			if(currentFrame>2){currentFrame=0;}
//		}
//		return score/Tools.max(1f, numKmers);
//	}
//	
//	/**
//	 * TODO
//	 * Evaluate the relative difference between left and right frequencies.
//	 * The purpose is to find locations where the left side looks noncoding and the right side looks coding.
//	 * Does not currently yield useful results.
//	 */
//	public static float scoreStart2(int point, byte[] bases, int stop, FrameStats stats){
//		final int k=stats.k;
//		
//		int start=point-45;
//		if(start<0 || stop>bases.length){return 0.5f;}
//
//		float left=calcKmerScore(start, Tools.min(point+k-2, bases.length), 0, bases, stats);
//		float right=calcKmerScore(point, stop-3, 0, bases, stats);
//		return right-left; //High numbers are likely to be starts; non-starts should be near 0.
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------           toString           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){

//		Collections.sort(fnames);
		taxIds.sort();
		
		bb.append("#BBMap "+Shared.BBTOOLS_VERSION_STRING+" Prokaryotic Gene Model\n");
		bb.append("#files");
		bb.tab().append(numFiles);
//		if(fnames.size()>5){
//			bb.tab().append(fnames.size());
//		}else{
//			for(String fname : fnames){
//				bb.tab().append(fname);
//			}
//		}
		bb.nl();
		bb.append("#taxIDs");
		for(int i=0; i<taxIds.size; i++){
			bb.tab().append(taxIds.get(i));
		}
		bb.nl();
//		bb.append("#k_inner\t").append(innerKmerLength).nl();
//		bb.append("#k_end\t").append(endKmerLength).nl();
//		bb.append("#start_left_offset\t").append(startLeftOffset).nl();
//		bb.append("#start_right_offset\t").append(startRightOffset).nl();
//		bb.append("#stop_left_offset\t").append(stopLeftOffset).nl();
//		bb.append("#stop_right_offset\t").append(stopRightOffset).nl();
		bb.append("#scaffolds\t").append(readsProcessed).nl();
		bb.append("#bases\t").append(basesProcessed).nl();
		bb.append("#genes\t").append(genesProcessed).nl();
		bb.append("#GC\t").append(gc(),2).nl();
		bb.append("#ACGTN");
		for(long x : baseCounts){
			bb.tab().append(x);
		}
		bb.nl();

		for(StatsContainer sc : allContainers){
			sc.appendTo(bb);
		}
		assert(allContainers.length>5) : allContainers.length;
		
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Stats             ----------------*/
	/*--------------------------------------------------------------*/

	public final StatsContainer statsCDS=new StatsContainer(CDS);
	public final StatsContainer statstRNA=new StatsContainer(tRNA);
	public final StatsContainer stats16S=new StatsContainer(r16S);
	public final StatsContainer stats23S=new StatsContainer(r23S);
	public final StatsContainer stats5S=new StatsContainer(r5S);
	public final StatsContainer stats18S=new StatsContainer(r18S);
	
	final StatsContainer[] rnaContainers=new StatsContainer[] {statstRNA, stats16S, stats23S, stats5S, stats18S};
	final StatsContainer[] allContainers=new StatsContainer[] {statsCDS, statstRNA, stats16S, stats23S, stats5S, stats18S};
	//public static int CDS=0, tRNA=1, r16S=2, r23S=3, r5S=4, r18S=5, r28S=6, RNA=7;
	
//	public final FrameStats innerKmerStats=new FrameStats("innerKmerStats", innerKmerLength, 3, 0);
//	public final FrameStats startStats=new FrameStats("startStats", endKmerLength, startFrames, startLeftOffset);
//	public final FrameStats stopStats=new FrameStats("stopStats", endKmerLength, stopFrames, stopLeftOffset);
//
//	public final FrameStats rrnaStartStats=new FrameStats("rrnaStart", 2, 16, 8);
//	public final FrameStats rrnaStopStats=new FrameStats("rrnaStop", 2, 16, 8);
//	
//	public final FrameStats trnaStats=new FrameStats("tRNA", rnaKmerLength, 1, 0);
//	public final FrameStats rrna16Sstats=new FrameStats("16S", rnaKmerLength, 1, 0);
//	public final FrameStats rrna23Sstats=new FrameStats("23S", rnaKmerLength, 1, 0);
//	public final FrameStats rrna5Sstats=new FrameStats("5S", rnaKmerLength, 1, 0);
//	public final FrameStats[] rnaKmerStats=new FrameStats[] {trnaStats, rrna16Sstats, rrna23Sstats, rrna5Sstats};
	
	/*--------------------------------------------------------------*/
	
//	public ArrayList<String> fnames=new ArrayList<String>();
	public int numFiles=0;
	public IntList taxIds=new IntList();

	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	long readsProcessed=0;
	long basesProcessed=0;
	long genesProcessed=0;
	long filesProcessed=0;
	long[] baseCounts=new long[5];
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Setters        ----------------*/
	/*--------------------------------------------------------------*/

	public synchronized void setStatics(){
//		assert(!setStatics);
		kInnerCDS=statsCDS.inner.k;
		kStartCDS=statsCDS.start.k;
		kStopCDS=statsCDS.stop.k;
		
		setStartLeftOffset(statsCDS.start.leftOffset);
		setStartRightOffset(statsCDS.start.rightOffset());
		
		setStopLeftOffset(statsCDS.stop.leftOffset);
		setStopRightOffset(statsCDS.stop.rightOffset());
		
		kInnerRNA=stats16S.inner.k;//TODO: Why is 16S used here?
		kStartRNA=stats16S.start.k;
		kStopRNA=stats16S.stop.k;
		setStatics=true;
	}
	
	public static void setInnerK(int k){
		kInnerCDS=k;
	}
	
	public static void setStartK(int k){
		kStartCDS=k;
	}
	
	public static void setStopK(int k){
		kStopCDS=k;
	}
	
	public static void setStartLeftOffset(int x){
		startLeftOffset=x;
		startFrames=startLeftOffset+startRightOffset+1;
//		System.err.println("startLeftOffset="+startLeftOffset+", startRightOffset="+startRightOffset+", frames="+startFrames);
	}
	
	public static void setStartRightOffset(int x){
		startRightOffset=x;
		startFrames=startLeftOffset+startRightOffset+1;
//		System.err.println("startLeftOffset="+startLeftOffset+", startRightOffset="+startRightOffset+", frames="+startFrames);
//		assert(false) : endLeftOffset+", "+endRightOffset+", "+endFrames;
	}
	
	public static void setStopLeftOffset(int x){
		stopLeftOffset=x;
		stopFrames=stopLeftOffset+stopRightOffset+1;
//		System.err.println("stopLeftOffset="+stopLeftOffset+", stopRightOffset="+stopRightOffset+", frames="+stopFrames);
	}
	
	public static void setStopRightOffset(int x){
		stopRightOffset=x;
		stopFrames=stopLeftOffset+stopRightOffset+1;
//		System.err.println("stopLeftOffset="+stopLeftOffset+", stopRightOffset="+stopRightOffset+", frames="+stopFrames);
//		assert(false) : endLeftOffset+", "+endRightOffset+", "+endFrames;
	}
	
	public static final boolean isStartCodon(int code){
		return code>=0 && code<=63 && isStartCodon[code];
	}
	public static final boolean isStopCodon(int code){
		return code>=0 && code<=63 && isStopCodon[code];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Class Init          ----------------*/
	/*--------------------------------------------------------------*/
	
	private static boolean[] makeIsCodon(String[] codons){
		boolean[] array=new boolean[64];
		for(String s : codons){
			int x=AminoAcid.toNumber(s);
			array[x]=true;
		}
		return array;
	}
	
	public static int kInnerCDS=6;
	public static int kStartCDS=3;
	public static int kStopCDS=3;
	
	static int startLeftOffset(){return startLeftOffset;}
	static int startRightOffset(){return startRightOffset;}
	static int startFrames(){return startFrames;}
	
	private static int startLeftOffset=21; //21 works well for k=4
	private static int startRightOffset=8; //10 works well for k=4
	private static int startFrames=startLeftOffset+startRightOffset+1;
	
	private static int stopLeftOffset=9;
	private static int stopRightOffset=12;
	private static int stopFrames=stopLeftOffset+stopRightOffset+1;
	
	private static boolean setStatics=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         More Statics         ----------------*/
	/*--------------------------------------------------------------*/
	
	//E. coli uses 83% AUG (3542/4284), 14% (612) GUG, 3% (103) UUG[7] and one or two others (e.g., an AUU and possibly a CUG).[8][9]
	public static String[] startCodons=new String[] {"ATG", "GTG", "TTG"};
	public static String[] extendedStartCodons=new String[] {"ATG", "GTG", "TTG", "ATT", "CTG", "ATA"};
	public static String[] stopCodons=new String[] {"TAG", "TAA", "TGA"};
	public static boolean[] isStartCodon=makeIsCodon(startCodons);
	public static boolean[] isStopCodon=makeIsCodon(stopCodons);
	
	/*--------------------------------------------------------------*/
	
	private static PrintStream outstream=System.err;
	public static boolean verbose=false;
	public static boolean errorState=false;
	
}
