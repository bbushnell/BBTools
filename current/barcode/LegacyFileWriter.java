package barcode;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map.Entry;
import java.util.Set;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import hiseq.IlluminaHeaderParser2;
import shared.LineParser1;
import shared.LineParserS1;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Writes legacy stats files when demultiplexing, for compatibility.
 * These are not very informative and not recommended for use.
 * @author BBushnell
 *
 */
public class LegacyFileWriter {

	public LegacyFileWriter(){}
	
//	Top_Unknown_Barcodes.csv
//	Lane,index,index2,# Reads,% of Unknown Barcodes,% of All Reads
//	7,AGGAGTTCTA,GGGGGGGGGG,3216744,0.017843,0.000839
//	7,TCGAATGATC,GGGGGGGGGG,2794420,0.015501,0.000729
//	7,CTTACCTGCT,GGGGGGGGGG,2709266,0.015028,0.000707
//	(1000 lines)
	long writeTopUnknownBarcodes(final Set<String> expectedSet, 
			HashMap<String, String> assignmentMap, final byte delimiter, 
			final String fname, final int lane, final int lines, 
			boolean overwrite) {
		if(fname==null) {return 0;}
		ArrayList<Barcode> list=new ArrayList<Barcode>(counts);
		Collections.sort(list);
		
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, false, true);
//		if(delimiter<=0) {
//			bsw.println("Lane,index,,# Reads,% of Unknown Barcodes,% of All Reads");
//		}else {
			bsw.println("Lane,index,index2,# Reads,% of Unknown Barcodes,% of All Reads");
//		}
		
		long total=0;
		long expected=0;
		for(Barcode bc : list) {
			total+=bc.count();
			if(assignmentMap.containsKey(bc.name)) {expected+=bc.count();}
		}
		final long unexpected=total-expected;
		final float invTotal=1f/Tools.max(1, total);
		final float invUnexpected=1f/Tools.max(1, unexpected);
		
//		System.err.println("total="+total+", expected="+expected+", unexpected="+unexpected);
		
		int printed=0;
		ByteBuilder bb=new ByteBuilder();
		for(Barcode bc : list) {
			if(!assignmentMap.containsKey(bc.name)) {
				if(printed>=lines) {break;}
				printed++;
				bb.clear().append(lane).comma();
				bc.appendIndex(bb, delimiter, 1).comma();
//				if(delimiter>0) {
					bc.appendIndex(bb, delimiter, 2).comma();//It's OK if there isn't one
//				}
//					System.err.println("count="+bc.count());
				bb.append(bc.count()).comma();
				bb.append(bc.count()*invUnexpected, 6).comma();
				bb.append(bc.count()*invTotal, 6).nl();
				bsw.print(bb);
			}
		}
		bsw.poison();
		return unexpected;
	}
	

////	Index_Hopping_Counts.csv
////	Lane,SampleID,index,index2,# Reads,% of Hopped Reads,% of All Reads
////	7,529965,CCGTGCGGTC,CCCGATATGG,55663469,,0.014523
////	7,529966,CACAGTTGTC,TAGGCAGGCT,52701428,,0.013750
////	7,529967,CCTCTTCAGA,GATCCAGTAG,48686653,,0.012703
////	...
////	7,,CCGTGCGGTC,TAGGCAGGCT,3821,0.000132,0.000001
////	7,,CCGTGCGGTC,GATCCAGTAG,5645,0.000196,0.000001
////	7,,CCGTGCGGTC,GCCTGGTAAG,4099,0.000142,0.000001
////	(all valid then all invalid pairs)
	
	
	void writeIndexHoppingCounts(final Set<String> expectedSet, 
			LinkedHashMap<String,String> sampleMap, final byte delimiter, final String fname, 
			final int lane, boolean overwrite) {
		if(delimiter<=0) {return;}
		LinkedHashSet<String> leftSet=new LinkedHashSet<String>();
		LinkedHashSet<String> rightSet=new LinkedHashSet<String>();
		LineParserS1 lp=new LineParserS1(delimiter);
		for(String s : expectedSet) {
			lp.set(s);
			assert(lp.terms()==2) : "Wrong number of barcodes in '"+s+
				"'; expected 2 for delimiter "+(int)delimiter;
			leftSet.add(lp.parseString(0));
			rightSet.add(lp.parseString(1));
		}
//		System.err.println(expectedSet.size());
		final HashMap<String, Barcode> countMap=new HashMap<String, Barcode>(Tools.max(8, counts.size()));
		float allMult=1f/Tools.max(1, (totalR1));
		for(Barcode bc : counts) {countMap.put(bc.name, bc);}
		
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, false, true);
		bsw.println("Lane,SampleID,index,index2,# Reads,% of Hopped Reads,% of All Reads");
		final ByteBuilder bb=new ByteBuilder();
		for(String s : expectedSet) {
			String sample=(sampleMap==null ? null : sampleMap.get(s));
			lp.set(s);
			bb.clear().append(lane).comma().append(sample==null ? "" : sample).comma();
			bb.append(lp.parseByteArray(0)).comma();
			bb.append(lp.parseByteArray(1)).comma();
			Barcode bc=countMap.get(s);
			long count=(bc==null ? 0 : bc.count());
			bb.append(count).comma().comma().append(count*allMult,6).nl();
			bsw.print(bb);
		}
//		bsw.poisonAndWait();
//		assert(false);
		long hoppedCount=0;
		for(String left : leftSet) {
			bb.clear().append(left).append(delimiter);
			final int len=bb.length;
			for(String right : rightSet) {
				bb.setLength(len).append(right);
				String pair=bb.toString();
				if(!expectedSet.contains(pair)) {
					Barcode bc=countMap.get(pair);
					if(bc!=null) {
						hoppedCount+=bc.count();
					}
				}
			}
		}
		final float hoppedMult=1f/(Tools.max(1, hoppedCount));
		
		ByteBuilder bb2=new ByteBuilder();
		for(String left : leftSet) {
			bb.clear().append(left).append(delimiter);
			final int len=bb.length;
			for(String right : rightSet) {
				bb.setLength(len).append(right);
				String pair=bb.toString();
				if(!expectedSet.contains(pair)) {
					Barcode bc=countMap.get(pair);
					long count=(bc==null ? 0 : bc.count());

					//String sample=(sampleMap==null ? null : sampleMap.get(s));
					lp.set(pair);
					bb2.clear().append(lane).comma().comma();
					bb2.append(lp.parseByteArray(0)).comma();
					bb2.append(lp.parseByteArray(1)).comma();
					bb2.append(count).comma().append(count*hoppedMult, 6).comma();
					bb2.append(count*allMult,6).nl();
					bsw.print(bb2);
				}
			}
		}
		bsw.poison();
	}
//		if(fname==null) {return 0;}
//		ArrayList<Barcode> list=new ArrayList<Barcode>(counts.size());
//		Collections.sort(list);
//		
//		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, false, true);
//		if(delimiter<=0) {
//			bsw.println("Lane,index,# Reads,% of Unknown Barcodes,% of All Reads");
//		}else {
//			bsw.println("Lane,index,index2,# Reads,% of Unknown Barcodes,% of All Reads");
//		}
//		
//		long total=0;
//		long expected=0;
//		for(Barcode bc : counts) {
//			total++;
//			if(expectedSet.contains(bc.name)) {expected+=bc.count();}
//		}
//		final long unexpected=total-expected;
//		final float invTotal=1f/Tools.max(1, total);
//		final float invUnexpected=1f/Tools.max(1, unexpected);
//		
//		int printed=0;
//		ByteBuilder bb=new ByteBuilder();
//		for(Barcode bc : counts) {
//			if(!expectedSet.contains(bc.name)) {
//				if(printed>=lines) {break;}
//				bb.clear().append(lane).comma();
//				bc.appendIndex(bb, delimiter, 1).comma();
//				if(delimiter>0) {bc.appendIndex(bb, delimiter, 2).comma();}
//				bb.append(bc.count()).comma();
//				bb.append(bc.count()*invUnexpected, 6).comma();
//				bb.append(bc.count()*invTotal, 6).nl();
//				bsw.print(bb);
//			}
//		}
//		bsw.poison();
//		return unexpected;
//	}
	
//	Quality_Metrics.csv
//	Lane,SampleID,index,index2,ReadNumber,Yield,YieldQ30,QualityScoreSum,Mean Quality Score (PF),% Q30
//	7,529965,CCGTGCGGTC,CCCGATATGG,1,8405183819,7964042904,327093323760,38.92,0.95
//	7,529965,CCGTGCGGTC,CCCGATATGG,2,8405183819,7794032223,323344091568,38.47,0.93
//	7,529966,CACAGTTGTC,TAGGCAGGCT,1,7957915628,7552613802,309911557188,38.94,0.95
//	...
//	7,Undetermined,,,1,27221786230,23126098298,1001525318468,36.79,0.85
//	7,Undetermined,,,2,27221786230,22727595251,990581252636,36.39,0.83
	int writeQualityMetrics(final Set<String> expectedSet, LinkedHashMap<String,String> sampleMap, 
			final byte delimiter, final String fname, final int lane, boolean overwrite) {
		if(fname==null) {return 0;}

		LinkedHashSet<String> set2=new LinkedHashSet<String>(expectedSet);
		set2.remove(UNDETERMINED);
		set2.add(UNDETERMINED);
		
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, false, true);
		bsw.println(QTag.header());
		ByteBuilder bb=new ByteBuilder();
		int lines=0;
		for(String name : set2) {
			String sample=(sampleMap==null ? null : sampleMap.get(name));
			QTag[] tags=tagMap.get(name);
			if(tags==null) {
				tags=new QTag[] {new QTag(), new QTag()};
			}
			for(int pairnum=0; pairnum<tags.length && tags[pairnum]!=null; pairnum++) {
				QTag tag=tags[pairnum];
				bb.clear();
				tag.appendTo(bb, lane, pairnum, sample==null ? UNDETERMINED : sample,
						delimiter, name==UNDETERMINED ? "" : name);
				bsw.print(bb.nl());
				lines++;
			}
		}
		bsw.poison();
		return lines;
	}
	
//	Demultiplex_Stats.csv
//	Lane,SampleID,Index,# Reads,# Perfect Index Reads,# One Mismatch Index Reads,# Two Mismatch Index Reads,% Reads,% Perfect Index Reads,% One Mismatch Index Reads,% Two Mismatch Index Reads
//	7,529965,CCGTGCGGTC-CCCGATATGG,55663469,49490086,6173383,0,0.0145,0.8891,0.1109,0.0000
//	7,529966,CACAGTTGTC-TAGGCAGGCT,52701428,48056108,4645320,0,0.0138,0.9119,0.0881,0.0000
//	7,529967,CCTCTTCAGA-GATCCAGTAG,48686653,44431620,4255033,0,0.0127,0.9126,0.0874,0.0000
	void writeDemultiplexStats(String fname, HashMap<String, String> assignmentMap,
			Set<String> expected, LinkedHashMap<String,String> sampleMap,
			int lane, boolean overwrite) {
		
		LinkedHashMap<String, long[]> hdistMap=new LinkedHashMap<String, long[]>(expected.size());
		final int maxHDist=3;
		for(String name : expected) {
			hdistMap.put(name, new long[maxHDist+1]);
		}
		
		HashMap<String, Barcode> countMap=new HashMap<String, Barcode>(counts.size());
		for(Barcode bc : counts) {countMap.put(bc.name, bc);}
		
		for(Entry<String, String> e : assignmentMap.entrySet()) {
			String key=e.getKey();
			String value=e.getValue();
			long[] distArray=hdistMap.get(value);
			assert(distArray!=null) : value;
			int dist=Tools.max(Barcode.hdistL(key, value), Barcode.hdistR(key, value));
			Barcode bc=countMap.get(key);
			long count=(bc==null ? 0 : bc.count());
			distArray[Tools.min(maxHDist, dist)]+=count;
		}
		
		String header="Lane,SampleID,Index,# Reads,# Perfect Index Reads,"
				+ "# One Mismatch Index Reads,# Two Mismatch Index Reads,"
				+ "% Reads,% Perfect Index Reads,% One Mismatch Index Reads,% Two Mismatch Index Reads";
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, false, true);
		bsw.println(header);
		ByteBuilder bb=new ByteBuilder();
		
		final float totalMult=1f/totalR1;
		for(String name : expected) {
			String sample=(sampleMap==null ? null : sampleMap.get(name));
			final long[] dist=hdistMap.get(name);
			final long sum=Tools.sum(dist);
			final float localMult=1f/sum;
			bb.clear().append(lane);
			bb.comma().append(sample==null ? "" : sample);
			bb.comma().append(name.replace('+', '-'));
			bb.comma().append(sum);
			for(int i=0; i<=2; i++) {bb.comma().append(dist[i]);}
			bb.comma().append(sum*totalMult, 4);
			for(int i=0; i<=2; i++) {bb.comma().append(dist[i]*localMult, 4);}
			bsw.print(bb.nl());
//			System.err.println(totalR1+", "+sum+", "+Arrays.toString(dist)+", "+(sum*totalMult)+", "+bb);
		}
		
		{
			long unknown=totalR1-totalR1Assigned;
			String sample="Undetermined";
//			final long[] dist=hdistMap.get(name);
			final long sum=unknown;
//			System.err.println(unknown+", "+totalR1+", "+totalR1Assigned);
			final long[] dist=new long[] {sum, 0, 0};
			final float localMult=1f/sum;
			bb.clear().append(lane);
			bb.comma().append(sample==null ? "" : sample);
			bb.comma().append("");
			bb.comma().append(sum);
			for(int i=0; i<=2; i++) {bb.comma().append(dist[i]);}
			bb.comma().append(sum*totalMult, 4);
			for(int i=0; i<=2; i++) {bb.comma().append(dist[i]*localMult, 4);}
			bsw.print(bb.nl());
		}
		
		bsw.poison();
	}
	
	
//	Demultiplex_Tile_Stats.csv
//	Lane,SampleID,Index,Tile,# Reads,% Reads
//	7,529965,CCGTGCGGTC-CCCGATATGG,1101,29931,0.0000
//	7,529965,CCGTGCGGTC-CCCGATATGG,1102,7526,0.0000
//	7,529965,CCGTGCGGTC-CCCGATATGG,1103,67763,0.0000


	void writeDemultiplexTileStats(String fname, LinkedHashMap<String,String> sampleMap, 
			Collection<String> expected, int lane, byte delimiter, boolean overwrite) {
		if(tileMap==null) {return;}
		if(expected==null) {expected=sampleMap.keySet();}
		String header="Lane,SampleID,Index,Tile,# Reads,% Reads";
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, false, true);
		bsw.println(header);
		ByteBuilder bb=new ByteBuilder();
		ArrayList<Integer> tiles=new ArrayList<Integer>(tileMap.keySet());
		Collections.sort(tiles);
		final float invTotal=1f/Tools.max(1, totalR1);
//		assert(false) : totalR1+", "+tileMap.size()+", "+expected.size();
		for(Integer tile : tiles) {
			LinkedHashMap<String, Barcode> tmap=tileMap.get(tile);
//			for(Entry<String, Barcode> e : tmap.entrySet()) {
//				Barcode bc=e.getValue();
			for(String name : expected) {
				Barcode bc=tmap.get(name);
				if(bc!=null) {
					String sampleID=(sampleMap==null ? null : sampleMap.get(bc.name));
					bb.clear();
					bb.append(lane);
					bb.comma().append(sampleID);
					Barcode.appendIndex(bb.comma(), delimiter, 1, bc.name);
					if(delimiter>=0 && (byte)bc.name.indexOf(delimiter)>0) {
						bb.dash();
						Barcode.appendIndex(bb, delimiter, 2, bc.name);
					}
					bb.comma().append(tile);
					bb.comma().append(bc.count());
					bb.comma().append(bc.count()*invTotal, 4).nl();
					bsw.print(bb);
				}
			}
		}
		bsw.poison();
	}
	
	
	void add(Read r, String barcode) {
		if(barcode==null) {barcode=UNDETERMINED;}
		QTag[] qtags=tagMap.get(barcode);
		if(qtags==null) {
			qtags=new QTag[] {new QTag(), new QTag()};
			tagMap.put(barcode, qtags);
		}
		qtags[r.pairnum()].add(r);
		String dest=(String) r.obj;
		if(r.pairnum()==0) {
			totalR1++;
			totalR1Assigned+=(barcode==UNDETERMINED ? 0 : 1);
			if(r.mate!=null) {
				totalR2++;
				totalR2Assigned+=(barcode==UNDETERMINED ? 0 : 1);
			}
//			assert(false) : (tileMap==null) +", "+dest+", "+barcode;
			if(tileMap!=null && dest!=null) {
				ihp.parse(r);
				Integer tile=ihp.tile();
				LinkedHashMap<String, Barcode> tmap=tileMap.get(tile);
				if(tmap==null) {
					tmap=new LinkedHashMap<String, Barcode>();
					tileMap.put(tile, tmap);
				}
				Barcode bc=tmap.get(dest);
				if(bc==null) {
					bc=new Barcode(dest);
					tmap.put(dest, bc);
				}
				assert(tileMap.size()>0 && tmap.size()>0) : 
					tileMap.size()+", "+tmap.size()+", "+tile+", "+dest+"\n"+ihp; 
				bc.increment(1);
			}
		}
	}
	

	LinkedHashMap<String, String> loadSampleMap(String fname, 
			byte barcodeDelimiter, boolean rcIndex1, boolean rcIndex2){
		if(fname==null) {return null;}
		LinkedHashMap<String, String> map=new LinkedHashMap<String, String>();
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		LineParser1 lp=new LineParser1(',');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			if(!Tools.startsWith(line, '#')) {
				lp.set(line);
				if(lp.terms()<2) {
					lp=new LineParser1('\t');
					lp.set(line);
					assert(lp.terms()>=2) : "Can't find delimiter in line '"+new String(line)+"'";
				}
				assert(lp.terms()>=2) : new String(line)+", "+lp;
				String key=lp.parseString(0), value=lp.parseString(1);
				key=BarcodeStats.fixBarcode(key, barcodeDelimiter, rcIndex1, rcIndex2);
				map.put(key, value);
			}
		}
		bf.close();
		return map;
	}
	
	static class QTag {
		
		int add(Read r) {
			return add(r.quality);
		}
		
		int add(byte[] quals) {
			reads++;
			bases+=quals.length;
			int q30=quals.length;
			int q30Count=0;
//			long sum=0;
			for(int i=0; i<quals.length; i++) {
				final int q=quals[i];
//				sum+=q;
				qsum.increment(i, q);
				qcount.increment(i, 1);
				if(q>=30) {
					q30Count++;
				}else if(i<q30){
					q30=i;
				}
			}
			q30PositionSum+=q30;
			q30CountSum+=q30Count;
			return q30;
		}
		
		long yieldBases() {return bases;}
		long yieldQ30() {return /*q30PositionSum*/q30CountSum;}
		long qualityScoreSum() {return (long)qsum.sum();}
		double meanQualityScore() {return qsum.sum()/bases;}
		double fractionQ30() {return yieldQ30()/(double)yieldBases();}
		
//		Lane,SampleID,index,index2,ReadNumber,Yield,YieldQ30,QualityScoreSum,Mean Quality Score (PF),% Q30
		ByteBuilder appendTo(ByteBuilder bb, int lane, int pairnum, 
				String sampleID, byte delimiter, String barcode) {
			bb.append(lane);
			bb.comma().append(sampleID);
			Barcode.appendIndex(bb.comma(), delimiter, 1, barcode);
			Barcode.appendIndex(bb.comma(), delimiter, 2, barcode);
			bb.comma().append(pairnum+1);
			bb.comma().append(yieldBases());
			bb.comma().append(yieldQ30());
			bb.comma().append(qualityScoreSum());
			bb.comma().append(meanQualityScore(),2);
			bb.comma().append(fractionQ30(),2);
			return bb;
		}
		
		static String header() {
			return "Lane,SampleID,index,index2,ReadNumber,"
					+ "Yield,YieldQ30,QualityScoreSum,Mean Quality Score (PF),% Q30";
		}
		
		long reads=0;
		long bases=0;
		long q30PositionSum=0;
		long q30CountSum=0;
		LongList qsum=new LongList();
		LongList qcount=new LongList();
	}

	private long totalR1=0;
	private long totalR2=0;
	private long totalR1Assigned=0;
	private long totalR2Assigned=0;
	
	Collection<Barcode> counts=null;
	final HashMap<String, QTag[]> tagMap=new HashMap<String, QTag[]>();
	private final static String UNDETERMINED="Undetermined";;
	LinkedHashMap<Integer, LinkedHashMap<String, Barcode>> tileMap=new LinkedHashMap<Integer, LinkedHashMap<String, Barcode>>();
	final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
	
	
}
