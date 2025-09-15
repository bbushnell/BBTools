package barcode;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import barcode.DemuxServer.Pair;
import fileIO.ByteStreamWriter;
import server.ServerTools;
import shared.KillSwitch;
import shared.LineParser2;
import shared.LineParserS1;
import shared.LineParserS2;
import shared.Shared;
import sketch.Sketch;
import structures.ByteBuilder;
import structures.StringNum;

/**
 * Sends data to the demux server.
 * @author BBushnell
 * @date June 24, 2024
 */
public class DemuxClient {
	
	public static void main(String[] args) {
		
		int coding=Sketch.RAW;
		if(args.length>0 && args[0].equalsIgnoreCase("A48")) {coding=Sketch.A48;}
		
		ArrayList<DemuxServer.Pair> pairs=new ArrayList<DemuxServer.Pair>();
		int len1, len2, delimiter;
		boolean dual=true;
		
		if(dual) {
			len1=10;
			len2=10;
			delimiter='+';
			pairs.add(new Pair("ACGTACGTAC+CGAACGAACG", "ACGTACGTAC+CGAACGAACG"));
			pairs.add(new Pair("ACGTACGTAC+NNNNNNNNNN", "ACGTACGTAC+CGAACGAACG"));
			pairs.add(new Pair("GGGGGGGGGG+NNNNNNNNNN", "GGGGGGGGGG+AGAACGAACG"));
			pairs.add(new Pair("GGGGGGGGGG+GNNNNNNNNN", "GGGGGGGGGG+CGAACGAACG"));
			pairs.add(new Pair("GGGGGGGGGG+TCGTACGTTT", "GGGGGGGGGG+CGAACGAACG"));
			pairs.add(new Pair("NNNNNNNNNN+AAAAAAAAAA", "GGGGGTGGGG+AAAAAAAAAA"));
		}else {

			len1=10;
			len2=0;
			delimiter=0;
			pairs.add(new Pair("ACGTACGTAC", "ACGTACGTAC"));
			pairs.add(new Pair("CCGTACGTAC", "ACGTACGTAC"));
			pairs.add(new Pair("AGGGGGGGGG", "GGGGGGGGGG"));
			pairs.add(new Pair("TGGGGGGGGG", "GGGGGGGGGG"));
			pairs.add(new Pair("GGGGGGGGGG", "GGGGGGGGGG"));
			pairs.add(new Pair("NNNNNNNNNN", "GGGGGTGGGG"));
		}
		
		Collections.sort(pairs);
		System.err.println("Expected:\n"+pairs);
		System.err.println("\nWriting:");
		
		ByteArrayOutputStream baos=new ByteArrayOutputStream();
		DemuxServer.writeToStream(pairs, baos, coding, len1, len2);
		try {
			baos.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String result=baos.toString();
		System.err.println("Created:\n"+result);
		

		System.err.println("Decoding.");
		HashMap<String, String> map=null;
		if(coding==Sketch.RAW) {
			map=toMapRaw(result, true);
		}else {
			map=toMapA48(result.getBytes(), true, coding, len1, len2, (byte)delimiter);
		}
		ArrayList<DemuxServer.Pair> pairs2=DemuxServer.toPairList(map);
		System.err.println("Decoded:\n"+pairs2);
	}
	
	public DemuxClient() {}
	
	public DemuxClient(String address_) {
		if(address_!=null) {address=address_;}
	}

	HashMap<String, String> getMap(DemuxData dd, boolean verbose) {
		String response=sendAndReceive(dd, verbose);
		if(dd.coding==Sketch.A48) {
			return toMapA48(response.getBytes(), verbose, dd);
		}else {
			return toMapRaw(response, verbose);
		}
	}
	
	private String sendAndReceive(DemuxData dd, boolean verbose) {
		if(verbose && false) {
			System.err.println("dd: type="+dd.type+", len1="+dd.length1+", len2="+dd.length2+
					", delimiter="+dd.barcodeDelimiter+", hdistsum="+dd.hdistSum+
					", expected="+dd.expectedList.size()+", counts="+dd.codeCounts.size());
		}
		
		ArrayList<byte[]> chunks = dd.encode();
		for(byte[] chunk : chunks) {
			assert(chunk[chunk.length-1]=='\n') : "'"+new String(chunk, chunk.length-60, 60)+"'";
		}

		if(verbose) {System.err.println("Sending "+chunks.size()+" chunks.");}

		if(writeSent!=null) {
			ByteStreamWriter bsw=ByteStreamWriter.makeBSW(writeSent, true, false, false);
			for(byte[] chunk : chunks) {
				bsw.print(chunk);
			}
			bsw.poisonAndWait();
		}
		
		StringNum response=ServerTools.sendAndReceive(chunks, address+"demux", verbose);

		//			System.err.println("Got response "+(result.s==null ? 0 : result.s.length()));
		if(!ServerTools.suppressErrors && (response.n<200 || response.n>299)){
			System.err.println("ERROR: Server returned code "+response.n+" and this message:\n"+response.s);
			KillSwitch.kill();
		}
		boolean e=checkForError(response.s);
		errorState|=e;
		if(e) {
			System.err.println("An error was encountered: \n"+response.s);
			KillSwitch.kill();
		}
		return response.s;
	}
	
	private static HashMap<String, String> toMapRaw(String response, boolean verbose){
		
		HashMap<String, String> map=new HashMap<String, String>();
		LineParserS2 responseParser=new LineParserS2('\n');
		LineParserS1 tp=new LineParserS1('\t');
		
		responseParser.set(response);
		while(responseParser.hasMore()) {
			responseParser.advance();
			if(responseParser.currentFieldLength()>0) {
				//System.err.println(responseParser.a()+", "+responseParser.b());
				String line=responseParser.parseStringFromCurrentField();
				tp.set(line);
				String key=null, value=null;
				try {
					key=tp.parseString(0);
					value=tp.parseString(1);
				} catch (Throwable t) {
					t.printStackTrace();
					System.err.println("'"+new String(line)+"'; "+tp.terms());
				}
				assert(!map.containsKey(key));
				map.put(key, value);
			}
		}
			
		return map;
	}
	
	private static HashMap<String, String> toMapA48(byte[] response, boolean verbose, DemuxData dd){
		return toMapA48(response, verbose, dd.coding, dd.length1, dd.length2, (byte)dd.barcodeDelimiter);
	}
	
	private static HashMap<String, String> toMapA48(byte[] response, boolean verbose, int coding, int len1, int len2, byte delimiter){
		HashMap<String, String> map=new HashMap<String, String>();
		LineParser2 responseParser=new LineParser2('\n');
		LineParser2 tp=new LineParser2('\t');
		final ByteBuilder bb=new ByteBuilder(32);
		
		byte[][] prev=new byte[4][];
		byte[][] current=new byte[4][];
		final int[] len=new int[] {len1, len2<1 ? len1 : len2, len1, len2};
		responseParser.set(response);
		while(responseParser.hasMore()) {
			responseParser.advance();
			if(responseParser.currentFieldLength()>0) {
				//System.err.println(responseParser.a()+", "+responseParser.b());
				byte[] line=responseParser.parseByteArrayFromCurrentField();
				tp.set(line);
				int term=0;
				while(tp.hasMore()) {
					tp.advance();
					if(tp.currentFieldLength()<1) {
						current[term]=prev[term];
					}else {
						current[term]=decodeA48(line, tp.a(), tp.b(), bb, len[term]);
					}
					term++;
				}
				
				assert(term==2 || term==4) : term+": '"+new String(line)+"'";
				final String key, value;
				if(term==2) {
					key=new String(current[0]);
					value=new String(current[1]);
				}else if(delimiter>0){
					key=bb.clear().append(current[0]).append(delimiter).append(current[1]).toString();
					value=bb.clear().append(current[2]).append(delimiter).append(current[3]).toString();
				}else {
					key=bb.clear().append(current[0]).append(current[1]).toString();
					value=bb.clear().append(current[2]).append(current[3]).toString();
				}
				assert(!map.containsKey(key));
				map.put(key, value);
			}
			byte[][] temp=prev;
			prev=current;
			current=temp;
		}
			
		return map;
	}
	
	public static byte[] decodeA48(byte[] term, int from, int to, ByteBuilder bb, int len) {
		bb.clear();
		long x=DemuxData.parseA48(term, from);
		DemuxData.decodeACGTN(x, bb, len);
		return bb.toBytes();
	}
	
	private static boolean checkForError(String s){
		if(s==null){return false;}
		return s.contains("java.io.IOException: Server returned HTTP response code:");
	}
	
	public boolean errorState=false;
	public String address=Shared.demuxServer();
	public static String writeSent=null;
	
}
