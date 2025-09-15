package barcode;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

import fileIO.ByteStreamWriter;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

public class BarcodeMappingStats {
	
	/*--------------------------------------------------------------*/
	/*----------------          Constructor         ----------------*/
	/*--------------------------------------------------------------*/
	
	public BarcodeMappingStats() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void merge(BarcodeMappingStats bs) {
		for(Entry<String, Barcode> e : bs.codeMap.entrySet()) {
			Barcode b=e.getValue();
			incrementCodeMap(b.name, b.count());
		}
		for(Entry<String, HashMap<String, Barcode>> e : bs.sourceMap.entrySet()) {
			final String readKey=e.getKey();
			final HashMap<String, Barcode> map=e.getValue();
			for(Entry<String, Barcode> ee : map.entrySet()) {
				final Barcode b=ee.getValue();
				final String refKey=ee.getKey();
				incrementSourceMap(readKey, refKey, b.count());
			}
		}
	}
	
	public void increment(Read r, String refKey){
		String barcode=r.barcode(true);
		incrementCodeMap(barcode, r.pairCount());
		incrementSourceMap(barcode, refKey==null ? "UNKNOWN" : refKey, r.pairCount());
	}
	
	public void incrementCodeMap(String key, long amt) {
		Barcode b=codeMap.get(key);
		if(b==null){
			b=new Barcode(key);
			codeMap.put(key, b);
		}
		b.increment(amt);
	}
	
	public void incrementSourceMap(String readKey, String refKey, long amt) {
		HashMap<String, Barcode> map=sourceMap.get(readKey);
		if(map==null){
			map=new HashMap<String, Barcode>();
			sourceMap.put(readKey, map);
		}
		Barcode b=map.get(refKey);
		if(b==null){
			b=new Barcode(refKey);
			map.put(refKey, b);
		}
		b.increment(amt);
	}

	public void writeStats(String outbarcodes, boolean overwrite) {
		ByteBuilder bb=new ByteBuilder();
		ArrayList<Barcode> codeList=toSortedList(codeMap);
		final long sum=sum(codeList);
		final double invSum=1.0/(Tools.max(1, sum));

		ByteStreamWriter bsw=new ByteStreamWriter(outbarcodes, overwrite, false, true);
		bsw.start();

		bsw.println("#Reads\t"+sum);
		bsw.println("#Barcode\tSource\tCount\tFraction");
		for(Barcode bc : codeList) {
			HashMap<String, Barcode> map=sourceMap.get(bc.name);
			ArrayList<Barcode> sourceList=toSortedList(map);
			final long sum2=sum(sourceList);
			final double invSum2=1.0/(Tools.max(1, sum2));
			for(Barcode source : sourceList) {
				bb.append(bc.name).tab().append(source.name).tab().append(source.count()).tab().append(source.count()*invSum2, 6).nl();
				bsw.print(bb);
				bb.clear();
			}
		}
		errorState|=bsw.poisonAndWait();
	}
	
	private static long sum(ArrayList<Barcode> list) {
		long sum=0;
		for(Barcode bc : list) {
			sum+=bc.count();
		}
		return sum;
	}
	
	private static ArrayList<Barcode> toSortedList(HashMap<String, Barcode> map){
		if(map==null || map.isEmpty()){return null;}
		ArrayList<Barcode> list=new ArrayList<Barcode>(map.size());
		for(Entry<String, Barcode> e : map.entrySet()) {
			list.add(e.getValue());
		}
		Collections.sort(list);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Raw counts of barcodes */
	public HashMap<String, Barcode> codeMap=new HashMap<String, Barcode>();
	
	/** 
	 * A table of tables.  Key1 is the barcode of a read; key2 is where the read mapped to.
	 * The barcode for key2 tracks the number of times reads with  key1 barcode mapped to key2's reference. 
	 * E.G. if key1 is ABC-DEF then the top barcode in its table would be expected to
	 * be ABC-DEF, and other entries would indicate contamination of that library.
	 * 
	 */
	public HashMap<String, HashMap<String, Barcode>> sourceMap=new HashMap<String, HashMap<String, Barcode>>();
	
	public boolean errorState=false;
	
}
