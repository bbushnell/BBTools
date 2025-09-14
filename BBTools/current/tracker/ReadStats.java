package tracker;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

import align2.QualityTools;
import barcode.Barcode;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import shared.Shared;
import shared.Tools;
import shared.TrimRead;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;
import structures.LongList;
import structures.SuperLongList;


/**
 * @author Brian Bushnell
 * @date Mar 18, 2013
 *
 */
public class ReadStats {
	
	public ReadStats(){this(true);}
		
	public ReadStats(boolean addToList){
		if(addToList){
			synchronized(ReadStats.class){
				objectList.add(this);
			}
		}

		if(COLLECT_QUALITY_STATS){
			aqualArray=new long[2][QMAX2];
			qualLength=new long[2][MAXLEN];
			qualSum=new long[2][MAXLEN];
			qualSumDouble=new double[2][MAXLEN];
			bqualHistOverall=new long[QMAX2];
		}else{
			aqualArray=null;
			qualLength=null;
			qualSum=null;
			qualSumDouble=null;
			bqualHistOverall=null;
		}
		
		if(BQUAL_HIST_FILE!=null || COLLECT_QUALITY_STATS || COLLECT_QUALITY_ACCURACY){
			bqualHist=new long[2][MAXLEN][QMAX2];
		}else{
			bqualHist=null;
		}
		
		if(QUAL_COUNT_HIST_FILE!=null){
			qcountHist=new long[2][QMAX2];
		}else{
			qcountHist=null;
		}

		if(COLLECT_MATCH_STATS){
			matchSum=new long[2][MAXLEN];
			delSum=new long[2][MAXLEN];
			insSum=new long[2][MAXLEN];
			subSum=new long[2][MAXLEN];
			nSum=new long[2][MAXLEN];
			clipSum=new long[2][MAXLEN];
			otherSum=new long[2][MAXLEN];
		}else{
			matchSum=null;
			delSum=null;
			insSum=null;
			subSum=null;
			nSum=null;
			clipSum=null;
			otherSum=null;
		}
		
		if(COLLECT_QUALITY_ACCURACY){
			qualMatch=new long[QMAX2];
			qualSub=new long[QMAX2];
			qualIns=new long[QMAX2];
			qualDel=new long[QMAX2];
		}else{
			qualMatch=null;
			qualSub=null;
			qualIns=null;
			qualDel=null;
		}
		
		if(COLLECT_INSERT_STATS){
			insertHist=new LongList(MAXLEN);
		}else{
			insertHist=null;
		}

		if(COLLECT_BASE_STATS){
			baseHist=new LongList[2][5];
			for(int i=0; i<baseHist.length; i++){
				for(int j=0; j<baseHist[i].length; j++){
					baseHist[i][j]=new LongList(400);
				}
			}
		}else{
			baseHist=null;
		}
		

		if(COLLECT_INDEL_STATS){
			insHist=new LongList(100);
			delHist=new LongList(100);
			delHist2=new LongList(100);
		}else{
			insHist=null;
			delHist=null;
			delHist2=null;
		}
		
		if(COLLECT_GC_STATS){
			gcHist=new long[GC_BINS+1];
		}else{
			gcHist=null;
		}
		
		if(COLLECT_ENTROPY_STATS){
			entropyHist=new long[ENTROPY_BINS+1];
			eTracker=new EntropyTracker(Shared.AMINO_IN, 0, true);
		}else{
			entropyHist=null;
			eTracker=null;
		}
		
		if(COLLECT_ERROR_STATS){
			errorHist=new LongList(100);
		}else{
			errorHist=null;
		}
		
		if(COLLECT_LENGTH_STATS){
			lengthHist=new SuperLongList(20000);
		}else{
			lengthHist=null;
		}
		
		if(COLLECT_IDENTITY_STATS){
			idHist=new long[ID_BINS+1];
			idBaseHist=new long[ID_BINS+1];
		}else{
			idHist=null;
			idBaseHist=null;
		}
		
		if(COLLECT_TIME_STATS){
			timeHist=new LongList(1001);
		}else{
			timeHist=null;
		}
		
		if(COLLECT_BARCODE_STATS) {
			barcodeMap=new HashMap<String, Barcode>();
		}else {
			barcodeMap=null;
		}
	}
	
	public static ReadStats mergeAll(){
		if(objectList==null || objectList.isEmpty()){return merged=null;}
		if(objectList.size()==1){return merged=objectList.get(0);}
		
		ReadStats x=new ReadStats(false);
		for(ReadStats rs : objectList){
			x.read2Count+=rs.read2Count;
			if(COLLECT_QUALITY_STATS){
				for(int i=0; i<MAXLEN; i++){
					x.qualLength[0][i]+=rs.qualLength[0][i];
					x.qualLength[1][i]+=rs.qualLength[1][i];
					x.qualSum[0][i]+=rs.qualSum[0][i];
					x.qualSum[1][i]+=rs.qualSum[1][i];
					x.qualSumDouble[0][i]+=rs.qualSumDouble[0][i];
					x.qualSumDouble[1][i]+=rs.qualSumDouble[1][i];
				}
				for(int i=0; i<x.aqualArray[0].length; i++){
					x.aqualArray[0][i]+=rs.aqualArray[0][i];
					x.aqualArray[1][i]+=rs.aqualArray[1][i];
				}
				for(int i=0; i<x.bqualHistOverall.length; i++){
					x.bqualHistOverall[i]+=rs.bqualHistOverall[i];
				}
				if(x.bqualHist!=null){
					for(int i=0; i<x.bqualHist.length; i++){
						for(int j=0; j<x.bqualHist[i].length; j++){
							for(int k=0; k<x.bqualHist[i][j].length; k++){
								x.bqualHist[i][j][k]+=rs.bqualHist[i][j][k];
							}
						}
					}
				}
				if(QUAL_COUNT_HIST_FILE!=null){
					for(int i=0; i<x.qcountHist.length; i++){
						for(int j=0; j<x.qcountHist[i].length; j++){
							x.qcountHist[i][j]+=rs.qcountHist[i][j];
						}
					}
				}
			}
			
			if(COLLECT_MATCH_STATS){
				for(int i=0; i<MAXLEN; i++){
					x.matchSum[0][i]+=rs.matchSum[0][i];
					x.matchSum[1][i]+=rs.matchSum[1][i];
					x.delSum[0][i]+=rs.delSum[0][i];
					x.delSum[1][i]+=rs.delSum[1][i];
					x.insSum[0][i]+=rs.insSum[0][i];
					x.insSum[1][i]+=rs.insSum[1][i];
					x.subSum[0][i]+=rs.subSum[0][i];
					x.subSum[1][i]+=rs.subSum[1][i];
					x.nSum[0][i]+=rs.nSum[0][i];
					x.nSum[1][i]+=rs.nSum[1][i];
					x.clipSum[0][i]+=rs.clipSum[0][i];
					x.clipSum[1][i]+=rs.clipSum[1][i];
					x.otherSum[0][i]+=rs.otherSum[0][i];
					x.otherSum[1][i]+=rs.otherSum[1][i];
				}
			}
			if(COLLECT_INSERT_STATS){
				x.insertHist.incrementBy(rs.insertHist);
				x.pairedCount+=rs.pairedCount;
				x.unpairedCount+=rs.unpairedCount;
			}
			if(COLLECT_BASE_STATS){
				for(int i=0; i<rs.baseHist.length; i++){
					for(int j=0; j<rs.baseHist[i].length; j++){
						x.baseHist[i][j].incrementBy(rs.baseHist[i][j]);
					}
				}
			}
			if(COLLECT_QUALITY_ACCURACY){
				for(int i=0; i<x.qualMatch.length; i++){
					x.qualMatch[i]+=rs.qualMatch[i];
					x.qualSub[i]+=rs.qualSub[i];
					x.qualIns[i]+=rs.qualIns[i];
					x.qualDel[i]+=rs.qualDel[i];
				}
			}
			

			if(COLLECT_INDEL_STATS){
				x.delHist.incrementBy(rs.delHist);
				x.delHist2.incrementBy(rs.delHist2);
				x.insHist.incrementBy(rs.insHist);
			}

			if(COLLECT_LENGTH_STATS){
				x.lengthHist.incrementBy(rs.lengthHist);
			}
			

			if(COLLECT_ERROR_STATS){
				x.errorHist.incrementBy(rs.errorHist);
			}
			
			if(COLLECT_GC_STATS){
				for(int i=0; i<rs.gcHist.length; i++){
					x.gcHist[i]+=rs.gcHist[i];
				}
			}
			
			if(COLLECT_ENTROPY_STATS){
				for(int i=0; i<rs.entropyHist.length; i++){
					x.entropyHist[i]+=rs.entropyHist[i];
				}
			}
			
			if(COLLECT_IDENTITY_STATS){
				for(int i=0; i<rs.idHist.length; i++){
					x.idHist[i]+=rs.idHist[i];
					x.idBaseHist[i]+=rs.idBaseHist[i];
				}
			}
			
			if(COLLECT_TIME_STATS){
				x.timeHist.incrementBy(rs.timeHist);
			}
			
			if(COLLECT_BARCODE_STATS){
				for(Entry<String, Barcode> e : rs.barcodeMap.entrySet()) {
					final String key=e.getKey();
					final Barcode brs=e.getValue();
					Barcode bx=x.barcodeMap.get(key);
					if(bx==null) {
						x.barcodeMap.put(key, brs);
					}else{
						bx.increment(brs);
					}
				}
			}

			x.gcMaxReadLen=Tools.max(x.gcMaxReadLen, rs.gcMaxReadLen);
			x.idMaxReadLen=Tools.max(x.idMaxReadLen, rs.idMaxReadLen);
		}
		
		merged=x;
		return x;
	}
	
	public void addToHistograms(Read r) {
		if(COLLECT_QUALITY_STATS){addToQualityHistogram(r);}
		if(COLLECT_BASE_STATS){addToBaseHistogram(r);}
		if(COLLECT_MATCH_STATS){addToMatchHistogram(r);}
		if(COLLECT_QUALITY_ACCURACY){addToQualityAccuracy(r);}
		
		if(COLLECT_ERROR_STATS){addToErrorHistogram(r);}
		if(COLLECT_INDEL_STATS){addToIndelHistogram(r);}
		if(COLLECT_LENGTH_STATS){addToLengthHistogram(r);}
		if(COLLECT_GC_STATS){addToGCHistogram(r);}
		if(COLLECT_ENTROPY_STATS){addToEntropyHistogram(r);}
		if(COLLECT_IDENTITY_STATS){addToIdentityHistogram(r);}
		if(COLLECT_BARCODE_STATS){addToBarcodeStats(r);}
	}
	
	public void addToQualityHistogram(final Read r){
		if(r==null){return;}
		addToQualityHistogram2(r);
		if(r.mate!=null){addToQualityHistogram2(r.mate);}
	}
	
	private void addToQualityHistogram2(final Read r){
		final int pairnum=r.samline==null ? r.pairnum() : r.samline.pairnum();
		if(r==null || r.quality==null || r.quality.length<1){return;}
		byte[] quals=r.quality, bases=r.bases;
		final Object obj=r.obj;
		if(obj!=null && obj.getClass()==TrimRead.class){
			quals=(pairnum==0 ? ((TrimRead)obj).qual1 : ((TrimRead)obj).qual2);
			bases=(pairnum==0 ? ((TrimRead)obj).bases1 : ((TrimRead)obj).bases2);
		}
		if(pairnum==1){read2Count++;}
		addToQualityHistogram(quals, pairnum);
		int x=Read.avgQualityByProbabilityInt(bases, quals, true, 0);
		aqualArray[pairnum][x]++;
		if(bqualHist!=null){
			addToBQualityHistogram(quals, pairnum);
		}
		if(QUAL_COUNT_HIST_FILE!=null){
			addToQCountHistogram(quals, pairnum);
		}
	}
	
	public void addToQualityHistogram(byte[] qual, int pairnum){
		if(qual==null || qual.length<1){return;}
		final int limit=Tools.min(qual.length, MAXLEN);
		final long[] ql=qualLength[pairnum], qs=qualSum[pairnum];
		final double[] qsd=qualSumDouble[pairnum];
		ql[limit-1]++;
		for(int i=0; i<limit; i++){
			qs[i]+=qual[i];
			qsd[i]+=QualityTools.PROB_ERROR[qual[i]];
		}
		for(byte q : qual){
			bqualHistOverall[q]++;
		}
	}
	
	private void addToBQualityHistogram(byte[] qual, int pairnum){
		if(qual==null || qual.length<1){return;}
		final int limit=Tools.min(qual.length, MAXLEN);
		final long[][] bqh=bqualHist[pairnum];
		for(int i=0; i<limit; i++){
			bqh[i][qual[i]]++;
		}
	}
	
	private void addToQCountHistogram(byte[] qual, int pairnum){
		if(qual==null || qual.length<1){return;}
		final long[] qch=qcountHist[pairnum];
		for(byte q : qual){
			qch[q]++;
		}
	}
	
	public void addToQualityAccuracy(final Read r){
		if(r==null){return;}
		addToQualityAccuracy(r, 0);
		if(r.mate!=null){addToQualityAccuracy(r.mate, 1);}
	}
	
	public void addToQualityAccuracy(final Read r, int pairnum){
		if(r==null || r.quality==null || r.quality.length<1 || !r.mapped() || r.match==null/* || r.discarded()*/){return;}
		final byte[] bases=r.bases;
		final byte[] qual=r.quality;
		byte[] match=r.match;
		
		if(r.shortmatch()){match=Read.toLongMatchString(match);}

		final boolean plus=(r.strand()==0);
		int bpos=0;
		byte lastm='A';
		for(int mpos=0; mpos<match.length/* && rpos<limit*/; mpos++){
			byte b=bases[bpos];
			byte q=qual[bpos];
			byte m=match[plus ? mpos : match.length-mpos-1];
			
			{
				if(m=='m'){
					qualMatch[q]++;
				}else if(m=='S'){
					qualSub[q]++;
				}else if(m=='I'){
					if(AminoAcid.isFullyDefined(b)){qualIns[q]++;}
				}else if(m=='N'){
					//do nothing
				}else if(m=='C'){
					//do nothing
				}else if(m=='V' || m=='i'){
					//do nothing
				}else if(m=='D'){
					if(lastm!=m){
						int x=bpos, y=bpos-1;
						if(x<qual.length){
							if(AminoAcid.isFullyDefined(bases[x])){
								qualDel[qual[x]]++;
							}
						}
						if(y>=0){
							if(AminoAcid.isFullyDefined(bases[y])){
								qualDel[qual[y]]++;
							}
						}
					}
					bpos--;
				}else if(m=='d'){
					bpos--;
				}else{
					assert(!Tools.isDigit(m)) : ((char)m);
				}
			}

			bpos++;
			lastm=m;
		}
		
	}
	
	public void addToErrorHistogram(final Read r){
		if(r==null){return;}
		addToErrorHistogram(r, 0);
		if(r.mate!=null){addToErrorHistogram(r.mate, 1);}
	}
	
	private void addToErrorHistogram(final Read r, int pairnum){
		if(r==null || r.bases==null || r.length()<1 || !r.mapped() || r.match==null/* || r.discarded()*/){return;}
		r.toLongMatchString(false);
		int x=r.countSubs();
		errorHist.increment(x, 1);
	}
	
	public void addToLengthHistogram(final Read r){
		if(r==null){return;}
		addToLengthHistogram(r, 0);
		if(r.mate!=null){addToLengthHistogram(r.mate, 1);}
	}
	
	private void addToLengthHistogram(final Read r, int pairnum){
		if(r==null || r.bases==null){return;}
		int x=r.length();//Tools.min(r.length(), MAXLENGTHLEN); Old style before SLL
		lengthHist.increment(x, 1);
	}
	
	public void addToGCHistogram(final Read r1){
		if(r1==null){return;}
		final Read r2=r1.mate;
		final int len1=r1.length(), len2=r1.mateLength();
		
		final float gc1=(len1>0 ? r1.gc() : -1);
		final float gc2=(len2>0 ? r2.gc() : -1);
		if(usePairGC){
			final float gc;
			if(r2==null){
				gc=gc1;
			}else{
				gc=(gc1*len1+gc2*len2)/(len1+len2);
			}
			addToGCHistogram(gc, len1+len2);
		}else{
			addToGCHistogram(gc1, len1);
			addToGCHistogram(gc2, len2);
		}
	}
	
	private void addToGCHistogram(final float gc, final int len){
		if(gc<0 || len<1){return;}
		gcHist[Tools.min(GC_BINS, (int)(gc*(GC_BINS+1)))]++;
		gcMaxReadLen=Tools.max(len, gcMaxReadLen);
	}
	
	public void addToEntropyHistogram(final Read r1){
		if(r1==null){return;}
		final Read r2=r1.mate;
		final int len1=r1.length(), len2=r1.mateLength();
		
		final float entropy1=(len1>0 ? eTracker.averageEntropy(r1.bases, allowEntropyNs) : -1);
		final float entropy2=(len2>0 ? eTracker.averageEntropy(r2.bases, allowEntropyNs) : -1);
		if(/* usePairEntropy */ false){
			final float entropy;
			if(r2==null){
				entropy=entropy1;
			}else{
				entropy=(entropy1*len1+entropy2*len2)/(len1+len2);
			}
			addToEntropyHistogram(entropy, len1+len2);
		}else{
			addToEntropyHistogram(entropy1, len1);
			addToEntropyHistogram(entropy2, len2);
		}
	}
	
	private void addToEntropyHistogram(final float entropy, final int len){
		if(entropy<0 || len<1){return;}
		entropyHist[Tools.min(ENTROPY_BINS, (int)(entropy*(ENTROPY_BINS+1)))]++;
	}
	
	public void addToIdentityHistogram(final Read r){
		if(r==null){return;}
		addToIdentityHistogram(r, 0);
		if(r.mate!=null){addToIdentityHistogram(r.mate, 1);}
	}
		
	private void addToIdentityHistogram(final Read r, int pairnum){
		if(r==null || r.bases==null || r.length()<1 || !r.mapped() || r.match==null/* || r.discarded()*/){return;}
		float id=r.identity();
		idHist[(int)(id*ID_BINS)]++;
		idBaseHist[(int)(id*ID_BINS)]+=r.length();
		idMaxReadLen=Tools.max(r.length(), idMaxReadLen);
	}
	
	public void addToTimeHistogram(final Read r){
		if(r==null){return;}
		addToTimeHistogram(r, 0);//Time for pairs is the same.
	}
	
	private void addToTimeHistogram(final Read r, int pairnum){
		if(r==null){return;}
		assert(r.obj!=null && r.obj.getClass()==Long.class);
		int x=(int)Tools.min(((Long)r.obj).longValue(), MAXTIMELEN);
		timeHist.increment(x, 1);
	}
	
	public boolean addToIndelHistogram(final Read r){
		if(r==null){return false;}
		boolean success=addToIndelHistogram(r, 0);
		if(r.mate!=null){success=addToIndelHistogram(r.mate, 1)|success;}
		return success;
	}

	/** Handles short match, long match, and reads with attached SamLines */
	private boolean addToIndelHistogram(final Read r, int pairnum){
		if(r==null || !r.mapped()){return false;}
		if(r.samline!=null){
			boolean success=addToIndelHistogram(r.samline);
			if(success){return true;}
		}
		if(r.match==null/* || r.discarded()*/){return false;}
		final byte[] match=r.match;

		byte lastLetter='?';
		boolean digit=false;
		int streak=0;
		for(int mpos=0; mpos<match.length; mpos++){
			final byte m=match[mpos];

			if(Tools.isDigit(m)){
				streak=streak*10+m-'0';
				digit=true;
			}else{
				if(lastLetter==m){
					streak++;
				}else{
					//						assert(streak>0 || (streak==0 && lastm=='?'));
					if(!digit){streak++;}
					digit=false;
					if(lastLetter=='D'){
						streak=Tools.min(streak, MAXDELLEN2);
						if(streak<MAXDELLEN){delHist.increment(streak, 1);}
						delHist2.increment(streak/100, 1);
						//						System.err.println("A. Del: "+streak+", "+MAXDELLEN+", "+MAXDELLEN2);
					}else if(lastLetter=='I'){
						streak=Tools.min(streak, MAXINSLEN);
						insHist.increment(streak, 1);
						//						System.err.println("A. Ins: "+streak+", "+MAXINSLEN);
					}
					streak=0;
				}
				lastLetter=m;
			}
		}
		
		{//Final symbol
			if(!digit){streak++;}
			digit=false;
			if(lastLetter=='D'){
				streak=Tools.min(streak, MAXDELLEN2);
				if(streak<MAXDELLEN){delHist.increment(streak, 1);}
				delHist2.increment(streak/100, 1);
				//						System.err.println("A. Del: "+streak+", "+MAXDELLEN+", "+MAXDELLEN2);
			}else if(lastLetter=='I'){
				streak=Tools.min(streak, MAXINSLEN);
				insHist.increment(streak, 1);
				//						System.err.println("A. Ins: "+streak+", "+MAXINSLEN);
			}
			streak=0;
		}
		return true;
	}

	private boolean addToIndelHistogram(SamLine sl){
		if(sl==null || sl.cigar==null || !sl.mapped()){
			return false;
		}
		final String cigar=sl.cigar;
//		final int pairnum=sl.pairnum();
		
		int count=0;
		for(int cpos=0; cpos<cigar.length(); cpos++){
			final char c=cigar.charAt(cpos);

			if(Tools.isDigit(c)){
				count=count*10+c-'0';
			}else{
				if(c=='I'){
					int streak=Tools.min(count, MAXINSLEN);
					insHist.increment(streak, 1);
//					System.err.println("A. Ins: "+streak+", "+MAXINSLEN);
				}else if(c=='D'){
					int streak=Tools.min(count, MAXDELLEN2);
					if(streak<MAXDELLEN){delHist.increment(streak, 1);}
					delHist2.increment(streak/100, 1);
//					System.err.println("A. Del: "+streak+", "+MAXDELLEN+", "+MAXDELLEN2);
				}else{
					//Ignore
				}
				count=0;
			}
		}
		assert(count==0) : count;
		return true;
	}
	
	public void addToMatchHistogram(final Read r){
		if(r==null){return;}
		addToMatchHistogram2(r);
		if(r.mate!=null){addToMatchHistogram2(r.mate);}
	}
	
	private void addToMatchHistogram2(final Read r){
		if(r==null || r.bases==null || r.length()<1 || !r.mapped() || r.match==null/* || r.discarded()*/){return;}
		final int pairnum=r.samline==null ? r.pairnum() : r.samline.pairnum();
		if(pairnum==1){read2Count++;}
		final byte[] bases=r.bases;
		final int limit=Tools.min(bases.length, MAXLEN);
		final long[] ms=matchSum[pairnum], ds=delSum[pairnum], is=insSum[pairnum],
				ss=subSum[pairnum], ns=nSum[pairnum], cs=clipSum[pairnum], os=otherSum[pairnum];
		
		byte[] match=r.match;
		if(r.shortmatch() && match!=null){match=Read.toLongMatchString(match);}
		
		if(match==null){
			for(int i=0; i<limit; i++){
				byte b=bases[i];
				if(b=='N'){ns[i]++;}
				else{os[i]++;}
			}
		}else{
			final boolean plus=(r.strand()==0);
			int bpos=0;
			byte lastm='A';
			for(int mpos=0; mpos<match.length && bpos<limit; mpos++){
				byte b=bases[bpos];//bases[plus ? rpos : bases.length-rpos-1];
				byte m=match[plus ? mpos : match.length-mpos-1];//match[mpos];
				if(b=='N'){
					if(m=='D'){
						if(lastm!=m){ds[bpos]++;}
						bpos--;
					}else{ns[bpos]++;}
				}else{
					if(m=='m'){
						ms[bpos]++;
					}else if(m=='S'){
						ss[bpos]++;
					}else if(m=='I'){
						is[bpos]++;
					}else if(m=='N' || m=='V'){
//						assert(false) : "\n"+r+"\n"+new String(Data.getChromosome(r.chrom).getBytes(r.start, r.stop))+"\nrpos="+rpos+", mpos="+mpos;
						os[bpos]++;
					}else if(m=='C'){
//						assert(false) : r;
						cs[bpos]++;
					}else if(m=='D'){
						if(lastm!=m){ds[bpos]++;}
						bpos--;
					}else if(m=='i'){
						os[bpos]++;
					}else if(m=='d'){
						bpos--;
					}else{
						os[bpos]++;
						assert(false) : "For read "+r.numericID+", unknown symbol in match string: ASCII "+m+" = "+(char)m;
					}
				}
				bpos++;
				lastm=m;
			}
		}
	}
	
	public void addToInsertHistogram(final Read r, boolean ignoreMappingStrand){
		if(verbose){
			System.err.print(r.numericID);
			if(r==null || r.mate==null || !r.mapped() || !r.mate.mapped() || !r.paired()){
				System.err.println("\n");
			}else{
				System.err.println("\t"+r.strand()+"\t"+r.insertSizeMapped(ignoreMappingStrand)+"\t"+r.mate.insertSizeMapped(ignoreMappingStrand));
			}
		}
		if(r==null || r.mate==null || !r.mapped() || !r.mate.mapped() || !r.paired()){
			unpairedCount++;
			return;
		}
		int x=Tools.min(MAXINSERTLEN, r.insertSizeMapped(ignoreMappingStrand));
		if(x>0){
			insertHist.increment(x, 1);
			pairedCount++;
		}else{
			unpairedCount++;
		}
//		assert(x!=1) : "\n"+r+"\n\n"+r.mate+"\n";
//		System.out.println("Incrementing "+x);
	}
	
	public void addToInsertHistogram(final SamLine r1){
		int x=r1.tlen;
		if(x<0) {x=-x;}
		x=Tools.min(MAXINSERTLEN, x);
		if(r1.pairedOnSameChrom() && x>0){
			pairedCount++;
			insertHist.increment(x, 1);
		}else {
			unpairedCount++;
		}
	}
	
	public void addToInsertHistogram(final SamLine r1, final SamLine r2){
		if(r1==null){return;}
		int x=insertSizeMapped(r1, r2, REQUIRE_PROPER_PAIR);
		if(verbose){
			System.err.println(r1.qname+"\t"+x);
		}
		x=Tools.min(MAXINSERTLEN, x);
		if(x>0){
			insertHist.increment(x, 1);
			pairedCount++;
		}else{
			unpairedCount++;
		}
	}
	
	/** This is untested and only gives approximate answers when overlapping reads contain indels.
	 * It may give incorrect answers for same-strange pairs that are shorter than read length.
	 * It might give negative answers but that would be a bug. */
	public static int insertSizeMapped(SamLine r1, SamLine r2, boolean requireProperPair){
		if(r2==null){return r1.length();}
		if(!r1.mapped() || !r2.mapped() || !r1.pairedOnSameChrom() || (requireProperPair && !r1.properPair())){
			return -1;
		}
		
		int a1=r1.start(true, false);
		int a2=r2.start(true, false);
		
		if(r1.strand()!=r2.strand()){
			if(r1.strand()==1){return insertSizeMapped(r2, r1, requireProperPair);}
		}else if(a1>a2){
			return insertSizeMapped(r2, r1, requireProperPair);
		}
		
		int b1=r1.stop(a1, true, false);
		int b2=r2.stop(a2, true, false);

		int clen1=r1.calcCigarLength(true, false);
		int clen2=r2.calcCigarLength(true, false);
		
		int mlen1=b1-a1+1;
		int mlen2=b2-a2+1;
		
		int dif1=mlen1-clen1;
		int dif2=mlen2-clen2;
		
		int mlen12=b2-a1+1;
		
		if(Tools.overlap(a1, b1, a2, b2)){//hard case
			return mlen12-Tools.max(dif1, dif2); //Approximate
		}else{//easy case
			return mlen12-dif1-dif2;
		}
	}
	
	public void addToBaseHistogram(final Read r){
		addToBaseHistogram2(r);
		if(r.mate!=null){addToBaseHistogram2(r.mate);}
	}
	
	public void addToBaseHistogram2(final Read r){
		if(r==null || r.bases==null){return;}
		final int pairnum=r.samline==null ? r.pairnum() : r.samline.pairnum();
		if(pairnum==1){read2Count++;}
		final byte[] bases=r.bases;
		final LongList[] lists=baseHist[pairnum];
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b]+1;
			lists[x].increment(i, 1);
		}
	}
	
	public void addToBarcodeStats(final Read r){
		String key=r.barcode(false);
		if(key==null) {key="NONE";}
		Barcode b=barcodeMap.get(key);
		if(b==null) {
			b=new Barcode(key);
			barcodeMap.put(key, b);
		}
		b.increment(r.pairCount());
//		assert(false) : barcodeMap;
	}
	
	public static boolean testFiles(boolean allowDuplicates){
		return Tools.testOutputFiles(overwrite, append, allowDuplicates,
				AVG_QUAL_HIST_FILE, QUAL_HIST_FILE, BQUAL_HIST_FILE, BQUAL_HIST_OVERALL_FILE, QUAL_COUNT_HIST_FILE,
				MATCH_HIST_FILE, INSERT_HIST_FILE, BASE_HIST_FILE, QUAL_ACCURACY_FILE, INDEL_HIST_FILE, ERROR_HIST_FILE, LENGTH_HIST_FILE,
				GC_HIST_FILE, ENTROPY_HIST_FILE, IDENTITY_HIST_FILE, TIME_HIST_FILE, BARCODE_STATS_FILE);
	}
	
	public static boolean writeAll(){
		if(collectingStats()){
			ReadStats rs=mergeAll();
			boolean paired=rs.read2Count>0;
			
			if(AVG_QUAL_HIST_FILE!=null){rs.writeAverageQualityToFile(AVG_QUAL_HIST_FILE, paired);}
			if(QUAL_HIST_FILE!=null){rs.writeQualityToFile(QUAL_HIST_FILE, paired);}
			if(BQUAL_HIST_FILE!=null){rs.writeBQualityToFile(BQUAL_HIST_FILE, paired);}
			if(BQUAL_HIST_OVERALL_FILE!=null){rs.writeBQualityOverallToFile(BQUAL_HIST_OVERALL_FILE);}
			if(QUAL_COUNT_HIST_FILE!=null){rs.writeQCountToFile(QUAL_COUNT_HIST_FILE, paired);}
			if(MATCH_HIST_FILE!=null){rs.writeMatchToFile(MATCH_HIST_FILE, paired);}
			if(INSERT_HIST_FILE!=null){rs.writeInsertToFile(INSERT_HIST_FILE);}
			if(BASE_HIST_FILE!=null){rs.writeBaseContentToFile(BASE_HIST_FILE, paired);}
			if(QUAL_ACCURACY_FILE!=null){rs.writeQualityAccuracyToFile(QUAL_ACCURACY_FILE);}
			
			if(INDEL_HIST_FILE!=null){rs.writeIndelToFile(INDEL_HIST_FILE);}
			if(ERROR_HIST_FILE!=null){rs.writeErrorToFile(ERROR_HIST_FILE);}
			if(LENGTH_HIST_FILE!=null){rs.writeLengthToFile(LENGTH_HIST_FILE);}
			if(GC_HIST_FILE!=null){rs.writeGCToFile(GC_HIST_FILE, true);}
			if(ENTROPY_HIST_FILE!=null){rs.writeEntropyToFile(ENTROPY_HIST_FILE, true);}
			if(IDENTITY_HIST_FILE!=null){rs.writeIdentityToFile(IDENTITY_HIST_FILE, true);}
			if(TIME_HIST_FILE!=null){rs.writeTimeToFile(TIME_HIST_FILE);}
			if(BARCODE_STATS_FILE!=null){rs.writeBarcodesToFile(BARCODE_STATS_FILE);}
			
			boolean b=rs.errorState;
			clear();
			return b;
		}
		clear();
		return false;
	}
	
	public void writeAverageQualityToFile(String fname, boolean writePaired){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		tsw.print("#Quality\tcount1\tfraction1"+(writePaired ? "\tcount2\tfraction2" : "")+"\n");

		long sum1=shared.Vector.sum(aqualArray[0]);
		long sum2=shared.Vector.sum(aqualArray[1]);
		double mult1=1.0/Tools.max(1, sum1);
		double mult2=1.0/Tools.max(1, sum2);
		
		long y=sum1+sum2;
		for(int i=0; i<aqualArray[0].length; i++){
			long x1=aqualArray[0][i];
			long x2=aqualArray[1][i];
			y-=x1;
			y-=x2;
			tsw.print(Tools.format("%d\t%d\t%.5f", i, x1, x1*mult1));
			if(writePaired){
				tsw.print(Tools.format("\t%d\t%.5f", x2, x2*mult2));
			}
			tsw.print("\n");
			if(y<=0){break;}
		}
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
	}
	
	public void writeQCountToFile(String fname, boolean writePaired){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		tsw.print("#Quality\tcount1\tfraction1"+(writePaired ? "\tcount2\tfraction2" : "")+"\n");

		long sum1=shared.Vector.sum(qcountHist[0]);
		long sum2=shared.Vector.sum(qcountHist[1]);
		double mult1=1.0/Tools.max(1, sum1);
		double mult2=1.0/Tools.max(1, sum2);
		
		long y=sum1+sum2;
		for(int i=0; i<qcountHist[0].length; i++){
			long x1=qcountHist[0][i];
			long x2=qcountHist[1][i];
			y-=x1;
			y-=x2;
			tsw.print(Tools.format("%d\t%d\t%.5f", i, x1, x1*mult1));
			if(writePaired){
				tsw.print(Tools.format("\t%d\t%.5f", x2, x2*mult2));
			}
			tsw.print("\n");
			if(y<=0){break;}
		}
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
	}
	
	public int q30(int pairnum) {
		return lengthAboveAverageQscore(30, pairnum);
	}
	
	public int lengthAboveAverageQscore(float limit, int pairnum) {
		final long[] qs=qualSum[pairnum];
		long[] ql=qualLength[pairnum];
		if(ql!=null) {ql=ql.clone();}//Prevents modification of originals
		
		for(int i=MAXLEN-2; i>=0; i--){ql[i]+=ql[i+1];}
		for(int i=0; i<MAXLEN; i++) {
			if(ql[i]<0) {return i;}
			float avg=qs[i]/(float)ql[i];
			if(avg<limit) {return i;}
		}
		return MAXLEN;
	}
	
	public void writeQualityToFile(String fname, boolean writePaired){
		StringBuilder sb=new StringBuilder();
		final boolean measure=(matchSum!=null);
		if(measure){
			if(writePaired){
				sb.append("#BaseNum\tRead1_linear\tRead1_log\tRead1_measured\tRead2_linear\tRead2_log\tRead2_measured\n");
			}else{
				sb.append("#BaseNum\tRead1_linear\tRead1_log\tRead1_measured\n");
			}
		}else{
			sb.append("#BaseNum\tRead1_linear\tRead1_log"+(writePaired ? "\tRead2_linear\tRead2_log" : "")+"\n");
		}
		
		final long[] qs1=qualSum[0], qs2=qualSum[1];
		long[] ql1=qualLength[0], ql2=qualLength[1];
		if(ql1!=null) {ql1=ql1.clone();}//Prevents modification of originals
		if(ql2!=null) {ql2=ql2.clone();}
		final double[] qsd1=qualSumDouble[0], qsd2=qualSumDouble[1];
		
		for(int i=MAXLEN-2; i>=0; i--){
			ql1[i]+=ql1[i+1];
			ql2[i]+=ql2[i+1];
		}

		double div1sum=0;
		double div2sum=0;
		double deviation1sum=0;
		double deviation2sum=0;
		
		if(writePaired){
			for(int i=0; i<MAXLEN && (ql1[i]>0 || ql2[i]>0); i++){
				final int a=i+1;
				double blin, clin, blog, clog;
				final double div1=(double)Tools.max(1, ql1[i]);
				final double div2=(double)Tools.max(1, ql2[i]);
				
				blin=qs1[i]/div1;
				clin=qs2[i]/div2;
				blog=qsd1[i]/div1;
				clog=qsd2[i]/div2;
				blog=QualityTools.probErrorToPhredDouble(blog);
				clog=QualityTools.probErrorToPhredDouble(clog);
				if(measure){
					double bcalc=calcQualityAtPosition(i, 0);
					double ccalc=calcQualityAtPosition(i, 1);

					div1sum+=div1;
					div2sum+=div2;
					deviation1sum+=Math.abs(blog-bcalc)*div1;
					deviation2sum+=Math.abs(clog-ccalc)*div2;
					
					sb.append(Tools.format("%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", a, blin, blog, bcalc, clin, clog, ccalc));
				}else{
					sb.append(Tools.format("%d\t%.3f\t%.3f\t%.3f\t%.3f\n", a, blin, blog, clin, clog));
				}
			}
		}else{
			for(int i=0; i<MAXLEN && ql1[i]>0; i++){
				final int a=i+1;
				double blin, blog;
				final double div1=(double)Tools.max(1, ql1[i]);
				
				blin=qs1[i]/div1;
				blog=qsd1[i]/div1;
				blog=QualityTools.probErrorToPhredDouble(blog);
				if(measure){
					double bcalc=calcQualityAtPosition(i, 0);

					div1sum+=div1;
					deviation1sum+=Math.abs(blog-bcalc)*div1;
					
					sb.append(Tools.format("%d\t%.3f\t%.3f\t%.3f\n", a, blin, blog, bcalc));
				}else{
					sb.append(Tools.format("%d\t%.3f\t%.3f\n", a, blin, blog));
				}
			}
		}

		final double stdev;
		final double variance;
		final double entropy;
		{
			double stdevSum=0;
			double varianceSum=0;
			double entropySum=0;
			long values=0;
			for(int pos=0; pos<MAXLEN; pos++) {
				long count=0;
				for(int rnum=0; rnum<2; rnum++) {
					final long[] array=bqualHist[rnum][pos];
					final long sum=Tools.sum(array);
					if(sum>0) {
						values++;;
						count+=sum;
						final double avg=Tools.averageHistogram(array);
						final double stdev0=Tools.standardDeviationHistogram(array);
						final double entropy0=calcEntropySuperSlow(array);
						entropySum+=entropy0;
						double varSum=0;
						for(int i=0; i<array.length; i++) {
							double var=Tools.absdif(i, avg);
							varSum+=(var*array[i]);
						}
						stdevSum+=stdev0;
						varianceSum+=(varSum/sum);
					}
				}
				if(count==0) {break;}
			}
			stdev=stdevSum/values;
			variance=varianceSum/values;
			entropy=entropySum/values;
		}
		
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		if(measure){
			if(writePaired){
				tsw.print(Tools.format("#Deviation\t%.4f\t%.4f\n", deviation1sum/div1sum, deviation2sum/div2sum));
			}else{
				tsw.print(Tools.format("#Deviation\t%.4f\n", deviation1sum/div1sum));
			}
			tsw.print(Tools.format("#Avg_STDev\t%.3f\n", stdev));
			tsw.print(Tools.format("#Diversity\t%.3f\n", variance));
			tsw.print(Tools.format("#Entropy\t%.3f\n", entropy));
		}
		tsw.print(sb);
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
	}
	
	private double calcQualityAtPosition(int pos, int pairnum){
		boolean includeNs=true;
		long m=matchSum[pairnum][pos];
		long d=delSum[pairnum][pos]; //left-adjacent deletion
		long d2=delSum[pairnum][Tools.min(pos, delSum[pairnum].length-1)]; //right-adjacent deletion
		long i=insSum[pairnum][pos];
		long s=subSum[pairnum][pos];
		long n=(includeNs ? nSum[pairnum][pos] : 0); //This only tracks no-calls, not no-refs.
		long good=Tools.max(0, m*2-d-d2+n/2);
		long total=Tools.max(0, m*2+i*2+s*2+n*2); //not d
		long bad=total-good;
		if(total<1){return 0;}
		double error=bad/(double)total;
		return QualityTools.probErrorToPhredDouble(error);
	}
	
	public void writeBQualityOverallToFile(String fname){
		final long[] cp30=Arrays.copyOf(bqualHistOverall, bqualHistOverall.length);
		for(int i=0; i<30; i++){cp30[i]=0;}
		
		final long sum=shared.Vector.sum(bqualHistOverall);
		final long median=Tools.percentileHistogram(bqualHistOverall, 0.5);
		final double mean=Tools.averageHistogram(bqualHistOverall);
		final double stdev=Tools.standardDeviationHistogram(bqualHistOverall);
		final double mean30=Tools.averageHistogram(cp30);
		final double stdev30=Tools.standardDeviationHistogram(cp30);
		final double mult=1.0/Tools.max(1, sum);
		long y=sum;
		
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		tsw.print("#Median\t"+median+"\n");
		tsw.print("#Mean\t"+Tools.format("%.3f", mean)+"\n");
		tsw.print("#STDev\t"+Tools.format("%.3f", stdev)+"\n");
		tsw.print("#Mean_30\t"+Tools.format("%.3f", mean30)+"\n");
		tsw.print("#STDev_30\t"+Tools.format("%.3f", stdev30)+"\n");
		tsw.print("#Quality\tbases\tfraction\n");
		
		for(int i=0; i<bqualHistOverall.length; i++){
			long x=bqualHistOverall[i];
			y-=x;
			tsw.print(Tools.format("%d\t%d\t%.5f", i, x, x*mult));
			tsw.print("\n");
			if(y<=0){break;}
		}
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
	}
	
	public void writeBQualityToFile(String fname, boolean writePaired){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		tsw.print("#BaseNum\tcount_1\tmin_1\tmax_1\tmean_1\tQ1_1\tmed_1\tQ3_1\tLW_1\tRW_1");
		if(writePaired){tsw.print("\tcount_2\tmin_2\tmax_2\tmean_2\tQ1_2\tmed_2\tQ3_2\tLW_2\tRW_2");}
		tsw.print("\n");

		for(int i=0; i<MAXLEN; i++){
			final long[] a1=bqualHist[0][i], a2=bqualHist[1][i];
			final long sum1=shared.Vector.sum(a1), sum2=shared.Vector.sum(a2);
			if(sum1<1 && sum2<1){break;}
			
			{
				final long a[]=a1, sum=sum1;
				
				final long weightedSum=Tools.sumHistogram(a);
				final long med=Tools.medianHistogram(a), min=Tools.minHistogram(a), max=Tools.maxHistogram(a);
				final long firstQuart=Tools.percentileHistogram(a, 0.25);
				final long thirdQuart=Tools.percentileHistogram(a, 0.75);
				final long leftWhisker=Tools.percentileHistogram(a, 0.02);
				final long rightWhisker=Tools.percentileHistogram(a, 0.98);
				final double mean=weightedSum*1.0/Tools.max(sum, 0);
				tsw.print(Tools.format("%d\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d", i, sum, min, max, mean, firstQuart, med, thirdQuart, leftWhisker, rightWhisker));
			}

			if(writePaired){
				final long a[]=a2, sum=sum2;
				
				final long weightedSum=Tools.sumHistogram(a);
				final long med=Tools.medianHistogram(a), min=Tools.minHistogram(a), max=Tools.maxHistogram(a);
				final long firstQuart=Tools.percentileHistogram(a, 0.25);
				final long thirdQuart=Tools.percentileHistogram(a, 0.75);
				final long leftWhisker=Tools.percentileHistogram(a, 0.02);
				final long rightWhisker=Tools.percentileHistogram(a, 0.98);
				final double mean=weightedSum*1.0/Tools.max(sum, 0);
				tsw.print(Tools.format("\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d", sum, min, max, mean, firstQuart, med, thirdQuart, leftWhisker, rightWhisker));
			}
			tsw.print("\n");
		}
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
	}
	
	public void writeQualityAccuracyToFile(String fname){
		
		int max=qualMatch.length;
		for(int i=max-1; i>=0; i--){
			if(qualMatch[i]+qualSub[i]+qualIns[i]+qualDel[i]>0){break;}
			max=i;
		}
		
		double devsum=0;
		double devsumSub=0;
		long observations=0;
		for(int i=0; i<max; i++){
			long qm=qualMatch[i]*2;
			long qs=qualSub[i]*2;
			long qi=qualIns[i]*2;
			long qd=qualDel[i];
			
			double phred=-1;
			double phredSub=-1;
			
			long sum=qm+qs+qi+qd;
			if(sum>0){
				double mult=1.0/sum;
				double subRate=(qs)*mult;
				double errorRate=(qs+qi+qd)*mult;
				
				phredSub=QualityTools.probErrorToPhredDouble(subRate);
				phred=QualityTools.probErrorToPhredDouble(errorRate);
				double deviation=phred-i;
				double deviationSub=phredSub-i;
				if(i==QMIN && deviation<0){deviation=0;}
				else if(i==QMAX && max==QMAX2 && deviation>0){deviation=0;}
				if(i==QMIN && deviationSub<0){deviationSub=0;}
				else if(i==QMAX && max==QMAX2 && deviationSub>0){deviationSub=0;}
				devsum+=(Math.abs(deviation)*sum);
				devsumSub+=(Math.abs(deviationSub)*sum);
				observations+=sum;
			}
		}
		
		double stdevSum=0;
		double varianceSum=0;
		double entropySum=0;
		long values=0;
		for(int pos=0; pos<MAXLEN; pos++) {
			long count=0;
			for(int rnum=0; rnum<2; rnum++) {
				final long[] array=bqualHist[rnum][pos];
				final long sum=Tools.sum(array);
				if(sum>0) {
					values++;;
					count+=sum;
					final double avg=Tools.averageHistogram(array);
					final double stdev=Tools.standardDeviationHistogram(array);
					final double entropy=calcEntropySuperSlow(array);
					entropySum+=entropy;
					double varSum=0;
					for(int i=0; i<array.length; i++) {
						double var=Tools.absdif(i, avg);
						varSum+=(var*array[i]);
					}
					stdevSum+=stdev;
					varianceSum+=(varSum/sum);
				}
			}
			if(count==0) {break;}
		}
		final double stdev=stdevSum/values;
		final double variance=varianceSum/values;
		final double entropy=entropySum/values;

		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		tsw.print(Tools.format("#Deviation\t%.3f\n", devsum/observations));
		tsw.print(Tools.format("#DeviationSub\t%.3f\n", devsumSub/observations));

		tsw.print(Tools.format("#Avg_STDev\t%.3f\n", stdev));
		tsw.print(Tools.format("#Diversity\t%.3f\n", variance));
		tsw.print(Tools.format("#Entropy\t%.3f\n", entropy));
		
		tsw.print("#Quality\tMatch\tSub\tIns\tDel\tTrueQuality\tTrueQualitySub\n");
		for(int i=0; i<max; i++){
			long qm=qualMatch[i]*2;
			long qs=qualSub[i]*2;
			long qi=qualIns[i]*2;
			long qd=qualDel[i];
			
			double phred=-1;
			double phredSub=-1;
			
			long sum=qm+qs+qi+qd;
			if(sum>0){
				double mult=1.0/sum;
				double subRate=(qs)*mult;
				double errorRate=(qs+qi+qd)*mult;
				
				phredSub=QualityTools.probErrorToPhredDouble(subRate);
				phred=QualityTools.probErrorToPhredDouble(errorRate);
				
//				System.err.println("sub: "+qs+"/"+sum+" -> "+subRate+" -> "+phredSub);
			}
			
			tsw.print(i+"\t"+qm+"\t"+qs+"\t"+qi+"\t"+qd);
			tsw.print(phred>=0 ? Tools.format("\t%.2f", phred) : "\t");
			tsw.print(phredSub>=0 ? Tools.format("\t%.2f\n", phredSub) : "\t\n");
			
//			System.err.println(qm+"\t"+qs+"\t"+qi+"\t"+qd);
		}
		
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
	}
	
	private float calcEntropySuperSlow(long[] counts){
		//Sum of entropy contributions from each term
		double eSum=0;
		
		final double total=Tools.sum(counts);
		
		//Loop over all nonzero kmer counts
		for(long count : counts){
			//Prevent NaN and INF
//			System.err.println("count="+count);
			if(count>0){
				//Fraction of total represented by this term
				double pk=count/(double)total;
//				System.err.println("pk="+pk);

				//Entropy contribution for this kmer's count
				double npklogpk=-pk*Math.log(pk);
//				System.err.println("npklogpk="+npklogpk);
				eSum+=npklogpk;
//				System.err.println("eSum="+eSum);
			}
		}
		//eSum now holds entropy in bits.
		
		//Multiplier to convert entropy to 0-1 scale.
		double multiplier=1/Math.log(total);
		
		//Adjust entropy to 0-1 scale based on window size
		float e=(float)(eSum*multiplier);
//		System.err.println("e="+e);
		assert(e>=0 && e<=1) : e+", "+eSum+", "+total+"\n"+Arrays.toString(counts)+"\n"+this;
		
		//Get rid of negative zero
		if(e<=0){e=0;}
		
		return e;
	}
	
	public void writeMatchToFile(String fname, boolean writePaired){
		if(!writePaired){
			writeMatchToFileUnpaired(fname);
			return;
		}
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#BaseNum\tMatch1\tSub1\tDel1\tIns1\tN1\tOther1\tMatch2\tSub2\tDel2\tIns2\tN2\tOther2\n");
		
		final long[] ms1=matchSum[0], ds1=delSum[0], is1=insSum[0],
				ss1=subSum[0], ns1=nSum[0], cs1=clipSum[0], os1=otherSum[0];
		final long[] ms2=matchSum[1], ds2=delSum[1], is2=insSum[1],
				ss2=subSum[1], ns2=nSum[1], cs2=clipSum[1], os2=otherSum[1];
		
		for(int i=0; i<MAXLEN; i++){
			int a=i+1;
			long sum1=ms1[i]+is1[i]+ss1[i]+ns1[i]+cs1[i]+os1[i]; //no deletions
			long sum2=ms2[i]+is2[i]+ss2[i]+ns2[i]+cs2[i]+os2[i]; //no deletions
			if(sum1==0 && sum2==0){break;}
			double inv1=1.0/(double)Tools.max(1, sum1);
			double inv2=1.0/(double)Tools.max(1, sum2);

			tsw.print(Tools.format("%d", a));
			tsw.print(Tools.format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f",
					ms1[i]*inv1, ss1[i]*inv1, ds1[i]*inv1, is1[i]*inv1, ns1[i]*inv1, (os1[i]+cs1[i])*inv1));
			tsw.print(Tools.format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f",
					ms2[i]*inv2, ss2[i]*inv2, ds2[i]*inv2, is2[i]*inv2, ns2[i]*inv2, (os2[i]+cs2[i])*inv2)
//					+", "+ms2[i]+", "+is2[i]+", "+ss2[i]+", "+ns2[i]+", "+cs2[i]+", "+os2[i]
					);
			tsw.print("\n");
		}
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
	}
	
	public void writeMatchToFileUnpaired(String fname){
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#BaseNum\tMatch1\tSub1\tDel1\tIns1\tN1\tOther1\n");
		
		final long[] ms1=matchSum[0], ds1=delSum[0], is1=insSum[0],
				ss1=subSum[0], ns1=nSum[0], cs1=clipSum[0], os1=otherSum[0];
		
		for(int i=0; i<MAXLEN; i++){
			int a=i+1;
			long sum1=ms1[i]+is1[i]+ss1[i]+ns1[i]+cs1[i]+os1[i]; //no deletions
			if(sum1==0){break;}
			double inv1=1.0/(double)Tools.max(1, sum1);

			tsw.print(Tools.format("%d", a));
			tsw.print(Tools.format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f",
					ms1[i]*inv1, ss1[i]*inv1, ds1[i]*inv1, is1[i]*inv1, ns1[i]*inv1, (os1[i]+cs1[i])*inv1)
//					+", "+ms1[i]+", "+is1[i]+", "+ss1[i]+", "+ns1[i]+", "+cs1[i]+", "+os1[i]
					);
			tsw.print("\n");
		}
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
	}
	
	public void writeInsertToFile(String fname){
		StringBuilder sb=new StringBuilder();
		sb.append("#Mean\t"+Tools.format("%.3f", Tools.averageHistogram(insertHist.array))+"\n");
		sb.append("#Median\t"+Tools.percentileHistogram(insertHist.array, 0.5)+"\n");
		sb.append("#Mode\t"+Tools.calcModeHistogram(insertHist.array)+"\n");
		sb.append("#STDev\t"+Tools.format("%.3f", Tools.standardDeviationHistogram(insertHist.array))+"\n");
		double percent=pairedCount*100.0/(pairedCount+unpairedCount);
		sb.append("#PercentOfPairs\t"+Tools.format("%.3f", percent)+"\n");
//		sb.append("#PercentOfPairs\t"+Tools.format("%.3f", matedPercent)+"\n");
		sb.append("#InsertSize\tCount\n");
		writeHistogramToFile(fname, sb.toString(), insertHist, !skipZeroInsertCount);
	}
	
	public void writeBaseContentToFile(String fname, boolean paired){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, false);
		bsw.start();
		bsw.print("#Pos\tA\tC\tG\tT\tN\n");
		
		int max=writeBaseContentToFile2(bsw, baseHist[0], 0);
		if(paired){
			writeBaseContentToFile2(bsw, baseHist[1], max);
		}
		
		bsw.poisonAndWait();
		errorState|=bsw.errorState;
	}
	
	private static int writeBaseContentToFile2(ByteStreamWriter bsw, LongList[] lists, int offset){
		int max=0;
		ByteBuilder sb=new ByteBuilder(100);
		for(LongList ll : lists){max=Tools.max(max, ll.size);}
		for(int i=0; i<max; i++){
			long a=lists[1].get(i);
			long c=lists[2].get(i);
			long g=lists[3].get(i);
			long t=lists[4].get(i);
			long n=lists[0].get(i);
			double mult=1.0/(a+c+g+t+n);

			sb.setLength(0);
			sb.append(i+offset).tab();
			sb.append(a*mult, 5).tab();
			sb.append(c*mult, 5).tab();
			sb.append(g*mult, 5).tab();
			sb.append(t*mult, 5).tab();
			sb.append(n*mult, 5);
			sb.nl();
			
			bsw.print(sb);
		}
		return max;
	}
	
	public void writeIndelToFile(String fname){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, false);
		bsw.start();
		bsw.print("#Length\tDeletions\tInsertions\n");
		
		int max=Tools.max(insHist.size, delHist.size);

		ByteBuilder bb=new ByteBuilder(100);
		for(int i=0; i<max; i++){
			long x=delHist.get(i);
			long y=insHist.get(i);
			if(x>0 || y>0 || !skipZeroIndel){
				bb.clear();
				bb.append(i).tab().append(x).tab().append(y).nl();
				bsw.print(bb);
			}
		}
		
		//TODO: Disabled because it was irritating when graphing.  Should write to a different file.
//		tsw.print("#Length_bin\tDeletions\n");
//		max=delHist2.size;
//		for(int i=0; i<max; i++){
//			long x=delHist2.get(i);
//			if(x>0 || !skipZeroIndel){
//				tsw.print((i*DEL_BIN)+"\t"+x+"\n");
//			}
//		}
		
		bsw.poisonAndWait();
		errorState|=bsw.errorState;
	}
	
	public void writeErrorToFile(String fname){
		writeHistogramToFile(fname, "#Errors\tCount\n", errorHist, false);
	}
	
	public void writeLengthToFile(String fname){
		writeHistogramToFile(fname, "#Length\tCount\n", lengthHist, false);
	}
	
	public void writeTimeToFile(String fname){
		writeHistogramToFile(fname, "#Time\tCount\n", timeHist, false);
	}
	
	public void writeBarcodesToFile(String fname) {
		ArrayList<Barcode> barcodes=new ArrayList<Barcode>();
		for(Barcode b : barcodeMap.values()) {
			barcodes.add(b);
		}
		Collections.sort(barcodes);
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, false);
		bsw.start();
		long sum=0;
		for(Barcode b : barcodes) {
			sum+=b.count();
		}
		bsw.print("#Reads\t").print(sum).nl();
		bsw.print("#Barcodes\t").print(barcodes.size()).nl();
		
		ByteBuilder bb=new ByteBuilder(40);
		for(Barcode b : barcodes) {
			bb.clear().append(b.name).tab().append(b.count()).nl();
			bsw.print(bb);
		}
		bsw.poisonAndWait();
		errorState|=bsw.errorState;
	}
	
	public void writeHistogramToFile(String fname, String header, LongList hist, boolean printZeros){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, false);
		bsw.start();
		bsw.print(header);
		
		int max=hist.size;

		ByteBuilder bb=new ByteBuilder(40);
		for(int i=0; i<max; i++){
			long x=hist.get(i);
			if(x>0 || printZeros){
				bb.clear().append(i).tab().append(x).nl();
				bsw.print(bb);
			}
		}
		bsw.poisonAndWait();
		errorState|=bsw.errorState;
	}
	
	public void writeHistogramToFile(String fname, String header, SuperLongList hist, boolean printZeros){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, false);
		bsw.start();
		bsw.print(header);
		
		hist.sort();
		long max=hist.max();
		long[] array=hist.array();
		LongList list=hist.list();

		ByteBuilder bb=new ByteBuilder(40);
		for(int i=0; i<array.length; i++){
			long x=array[i];
			if(x>0 || printZeros){
				bb.clear().append(i).tab().append(x).nl();
				bsw.print(bb);
			}
		}
		{//TODO: This part does not bother printing zeros regardless of the flag.
			long prev=-1;
			int streak=0;
			for(int i=0; i<list.size(); i++){
				long x=list.get(i);
				if(x==prev){streak++;}
				else{
					if(streak>0){//Print previous key
						bb.clear().append(prev).tab().append(streak).nl();
						bsw.print(bb);
					}
					streak=1;
					prev=x;
				}
			}
			if(streak>0){//Print final key
				bb.clear().append(prev).tab().append(streak).nl();
				bsw.print(bb);
			}
		}
		bsw.poisonAndWait();
		errorState|=bsw.errorState;
	}
	
	public void writeGCToFile(String fname, boolean printZeros){
		final long[] hist;
		if(GC_BINS_AUTO && gcMaxReadLen+1<gcHist.length){
			hist=Tools.downsample(gcHist, gcMaxReadLen+1);
		}else{
			hist=gcHist;
		}
		final int bins=hist.length;
		final double gcMult=100.0/Tools.max(1, bins-1);
		final long total=shared.Vector.sum(hist);
		final long max=Tools.max(hist);
		final double countsPerX=Tools.max(1, ((max*1000.0)/40));
		final double fractionMult=1.0/Tools.max(1, total);
		long sum=0;	
		
		GCMean=Tools.averageHistogram(hist)*gcMult;
		GCMedian=Tools.percentileHistogram(hist, 0.5)*gcMult;
		GCMode=Tools.calcModeHistogram(hist)*gcMult;
		GCSTDev=Tools.standardDeviationHistogram(hist)*gcMult;
		
		ByteBuilder bb=new ByteBuilder(256);
		bb.append("#Mean\t").append(GCMean, 3).nl();
		bb.append("#Median\t").append(GCMedian, 3).nl();
		bb.append("#Mode\t").append(GCMode, 3).nl();
		bb.append("#STDev\t").append(GCSTDev, 3).nl();
		if(GC_PLOT_X){
			bb.append("#GC\tCount\tCumulative\tPlot\n");
		}else{
			bb.append("#GC\tCount\n");
		}
		

//		bsw.print("#Mean\t"+Tools.format("%.3f", GCMean)+"\n");
//		bsw.print("#Median\t"+Tools.format("%.3f", GCMedian)+"\n");
//		bsw.print("#Mode\t"+Tools.format("%.3f", GCMode)+"\n");
//		bsw.print("#STDev\t"+Tools.format("%.3f", GCSTDev)+"\n");
//		if(GC_PLOT_X){
//			bsw.print("#GC\tCount\tCumulative\tPlot\n");
//		}else{
//			bsw.print("#GC\tCount\n");
//		}
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, false);
		bsw.start();
		bsw.print(bb);
		
		for(int i=0; i<bins; i++){
			long x=hist[i];
			sum+=x;
			if(x>0 || printZeros){
				//This assumes GC_BINS==100
//				tsw.print(i+"\t"+x+"\n");
				if(GC_PLOT_X){
					bb.clear();
					bb.append(i*gcMult, 1).tab().append(x).tab();
					bb.append(sum*fractionMult, 3).tab();
					
					int len=(int)((x*1000)/countsPerX);
					for(int j=0; j<len; j++){bb.append('X');}
					if(len<1 && x>0){
						if((x*1000f)/countsPerX>0.1f){bb.append('x');}
						else{bb.append('.');}
					}
					
					bb.append('\n');
					bsw.print(bb);
				}else{
					bb.clear().append(i*gcMult, 1).tab().append(x).nl();
					bsw.print(bb);
				}
			}
		}
		bsw.poisonAndWait();
		errorState|=bsw.errorState;
	}
	
	public void writeEntropyToFile(String fname, boolean printZeros){
		final long[] hist=entropyHist;
		
		final int bins=hist.length;
		final double mult=1.0/Tools.max(1, bins-1);
		final long total=shared.Vector.sum(hist);
		final long max=Tools.max(hist);
		final double countsPerX=Tools.max(1, ((max*1000.0)/40));
		final double fractionMult=1.0/Tools.max(1, total);
		long sum=0;	
		
		double mean=Tools.averageHistogram(hist)*mult;
		double median=Tools.percentileHistogram(hist, 0.5)*mult;
		double mode=Tools.calcModeHistogram(hist)*mult;
		double stdev=Tools.standardDeviationHistogram(hist)*mult;
		
		ByteBuilder bb=new ByteBuilder(256);
		bb.append("#Mean\t").append(mean, 6).nl();
		bb.append("#Median\t").append(median, 6).nl();
		bb.append("#Mode\t").append(mode, 6).nl();
		bb.append("#STDev\t").append(stdev, 6).nl();
		bb.append("#Value\tCount\n");
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, false);
		bsw.start();
		bsw.print(bb);
		
		for(int i=0; i<bins; i++){
			long x=hist[i];
			sum+=x;
			if(x>0 || printZeros){
				bb.clear().append(i*mult, 4).tab().append(x).nl();
				bsw.print(bb);
			}
		}
		bsw.poisonAndWait();
		errorState|=bsw.errorState;
	}
	
	public void writeIdentityToFile(String fname, boolean printZeros){
		final long[] hist, histb;
		if(ID_BINS_AUTO && idMaxReadLen+1<idHist.length){
			hist=Tools.downsample(idHist, idMaxReadLen+1);
			histb=Tools.downsample(idBaseHist, idMaxReadLen+1);
		}else{
			hist=idHist;
			histb=idBaseHist;
		}
		final int max=hist.length;
		final double mult=100.0/(max-1);
		
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		
		
		tsw.print("#Mean_reads\t"+Tools.format("%.3f", (Tools.averageHistogram(hist)*mult))+"\n");
		tsw.print("#Mean_bases\t"+(Tools.format("%.3f", Tools.averageHistogram(histb)*mult))+"\n");
		tsw.print("#Median_reads\t"+(int)Math.round(Tools.percentileHistogram(hist, 0.5)*mult)+"\n");
		tsw.print("#Median_bases\t"+(int)Math.round(Tools.percentileHistogram(histb, 0.5)*mult)+"\n");
		tsw.print("#Mode_reads\t"+(int)Math.round(Tools.calcModeHistogram(hist)*mult)+"\n");
		tsw.print("#Mode_bases\t"+(int)Math.round(Tools.calcModeHistogram(histb)*mult)+"\n");
		tsw.print("#STDev_reads\t"+Tools.format("%.3f", (Tools.standardDeviationHistogram(hist)*mult))+"\n");
		tsw.print("#STDev_bases\t"+Tools.format("%.3f", (Tools.standardDeviationHistogram(histb)*mult))+"\n");
		tsw.print("#Identity\tReads\tBases\n");
		
		for(int i=0; i<max; i++){
			long x=hist[i], x2=histb[i];
			if(x>0 || printZeros){
				tsw.print(Tools.format("%.1f", i*mult)+"\t"+x+"\t"+x2+"\n");
			}
		}
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
	}
	
	//Tracks to see if read2s have been encountered, for displaying stats.
	private long read2Count=0;

	public long pairedCount=0;
	public long unpairedCount=0;
	
	private final int QMIN=Read.MIN_CALLED_QUALITY();
	private final int QMAX=Read.MAX_CALLED_QUALITY();
	private final int QMAX2=QMAX+1;
	
	public final long[][] aqualArray;
	public final long[][] qualLength;
	public final long[][] qualSum;
	
	public final long[][][] bqualHist;
	public final long[] bqualHistOverall;
	
	public final long[][] qcountHist;
	
	public final double[][] qualSumDouble;
	
	public final long[][] matchSum;
	public final long[][] delSum;
	public final long[][] insSum;
	public final long[][] subSum;
	public final long[][] nSum;
	public final long[][] clipSum;
	public final long[][] otherSum;

	public final long[] qualMatch;
	public final long[] qualSub;
	public final long[] qualIns;
	public final long[] qualDel;

	public final long[] gcHist;
	public final long[] entropyHist;
	public final EntropyTracker eTracker;
	public final long[] idHist;
	public final long[] idBaseHist;
	private int gcMaxReadLen=1;
	private int idMaxReadLen=1;

	public final LongList[][] baseHist;
	
	public final HashMap<String, Barcode> barcodeMap;
	
	/** Insert size */
	public final LongList insertHist;
	/** Read length */
	public final SuperLongList lengthHist;
	/** Number errors per read */
	public final LongList errorHist;
	/** Insertion length */
	public final LongList insHist;
	/** Deletion length */
	public final LongList delHist;
	/** Deletion length, binned  */
	public final LongList delHist2;
	/** Time */
	public final LongList timeHist;
	
	public static boolean REQUIRE_PROPER_PAIR=true;
	public static int MAXLEN=6000;
	public static int MAXINSERTLEN=80000;
	@Deprecated
	public static int MAXLENGTHLEN=80000;
	public static final int MAXTIMELEN=80000;
	public static final int MAXINSLEN=1000;
	public static final int MAXDELLEN=1000;
	public static final int MAXDELLEN2=1000000;
	public static final int DEL_BIN=100;
	public static int ID_BINS=100;
	public static boolean ID_BINS_AUTO=false;
	public static int GC_BINS=100;
	public static boolean GC_BINS_AUTO=false;
	public static boolean GC_PLOT_X=false;
	public static int ENTROPY_BINS=1000;

	public static double GCMean;
	public static double GCMedian;
	public static double GCMode;
	public static double GCSTDev;

//	public static double entropyMean;
//	public static double entropyMedian;
//	public static double entropyMode;
//	public static double entropySTDev;
	
	public boolean errorState=false;
	
	public static ReadStats merged=null;
	
//	public static double matedPercent=0;
	
	public static void clear(){
		COLLECT_QUALITY_STATS=false;
		COLLECT_QUALITY_ACCURACY=false;
		COLLECT_MATCH_STATS=false;
		COLLECT_INSERT_STATS=false;
		COLLECT_BASE_STATS=false;
		COLLECT_INDEL_STATS=false;
		COLLECT_GC_STATS=false;
		COLLECT_ENTROPY_STATS=false;
		COLLECT_ERROR_STATS=false;
		COLLECT_LENGTH_STATS=false;
		COLLECT_IDENTITY_STATS=false;
		COLLECT_TIME_STATS=false;
		COLLECT_BARCODE_STATS=false;
		
		AVG_QUAL_HIST_FILE=null;
		QUAL_HIST_FILE=null;
		BQUAL_HIST_FILE=null;
		QUAL_COUNT_HIST_FILE=null;
		BQUAL_HIST_OVERALL_FILE=null;
		QUAL_ACCURACY_FILE=null;
		MATCH_HIST_FILE=null;
		INSERT_HIST_FILE=null;
		BASE_HIST_FILE=null;
		INDEL_HIST_FILE=null;
		ERROR_HIST_FILE=null;
		LENGTH_HIST_FILE=null;
		GC_HIST_FILE=null;
		ENTROPY_HIST_FILE=null;
		IDENTITY_HIST_FILE=null;
		TIME_HIST_FILE=null;
		BARCODE_STATS_FILE=null;
		if(objectList!=null) {objectList.clear();}
//		objectList=null;
	}
	
	public static ArrayList<ReadStats> objectList=new ArrayList<ReadStats>();
	public static boolean COLLECT_QUALITY_STATS=false;
	public static boolean COLLECT_QUALITY_ACCURACY=false;
	public static boolean COLLECT_MATCH_STATS=false;
	public static boolean COLLECT_INSERT_STATS=false;
	public static boolean COLLECT_BASE_STATS=false;
	public static boolean COLLECT_INDEL_STATS=false;
	public static boolean COLLECT_GC_STATS=false;
	public static boolean COLLECT_ENTROPY_STATS=false;
	public static boolean COLLECT_ERROR_STATS=false;
	public static boolean COLLECT_LENGTH_STATS=false;
	public static boolean COLLECT_IDENTITY_STATS=false;
	public static boolean COLLECT_TIME_STATS=false;
	public static boolean COLLECT_BARCODE_STATS=false;
	
	public static boolean collectingStats(){
		return COLLECT_QUALITY_STATS || COLLECT_QUALITY_ACCURACY || COLLECT_MATCH_STATS || COLLECT_INSERT_STATS
				|| COLLECT_BASE_STATS || COLLECT_INDEL_STATS || COLLECT_GC_STATS || COLLECT_ENTROPY_STATS
				|| COLLECT_ERROR_STATS || COLLECT_LENGTH_STATS || COLLECT_IDENTITY_STATS || COLLECT_TIME_STATS
				|| COLLECT_BARCODE_STATS;
	}
	
	public static boolean usePairGC=true;
	public static boolean allowEntropyNs=true;
	
	public static String AVG_QUAL_HIST_FILE=null;
	public static String QUAL_HIST_FILE=null;
	public static String BQUAL_HIST_FILE=null;
	public static String QUAL_COUNT_HIST_FILE=null;
	public static String BQUAL_HIST_OVERALL_FILE=null;
	public static String QUAL_ACCURACY_FILE=null;
	public static String MATCH_HIST_FILE=null;
	public static String INSERT_HIST_FILE=null;
	public static String BASE_HIST_FILE=null;
	public static String INDEL_HIST_FILE=null;
	public static String ERROR_HIST_FILE=null;
	public static String LENGTH_HIST_FILE=null;
	public static String GC_HIST_FILE=null;
	public static String ENTROPY_HIST_FILE=null;
	public static String IDENTITY_HIST_FILE=null;
	public static String TIME_HIST_FILE=null;
	public static String BARCODE_STATS_FILE=null;
	
	public static boolean overwrite=true;
	public static boolean append=false;
	public static final boolean verbose=false;

	public static boolean skipZeroInsertCount=true;
	public static boolean skipZeroIndel=true;
	
}
