package barcode;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import barcode.DemuxServer.Pair;
import fileIO.ByteStreamWriter;
import parse.LineParser2;
import parse.LineParserS1;
import parse.LineParserS2;
import server.ServerTools;
import shared.KillSwitch;
import shared.Shared;
import sketch.Sketch;
import structures.ByteBuilder;
import structures.StringNum;

/**
 * Sends barcode-count requests to a demux server and decodes correction maps
 * returned by the server.
 * <p>
 * Requests are serialized by {@link DemuxData}; responses are either raw
 * tab-delimited pairs or A48-compressed pairs produced by
 * {@link DemuxServer#writeToStream(ArrayList, java.io.OutputStream, int, int, int)}.
 *
 * @author BBushnell
 * @date June 24, 2024
 */
public class DemuxClient {

/**
	 * Diagnostic harness for exercising the server response encoders and the
	 * matching local decoders with sample barcode pairs.
	 * @param args Optional encoding type ("A48" for compressed encoding)
	 */
	public static void main(String[] args){

		int coding=Sketch.RAW;
		if(args.length>0 && args[0].equalsIgnoreCase("A48")){coding=Sketch.A48;}

		ArrayList<DemuxServer.Pair> pairs=new ArrayList<DemuxServer.Pair>();
		int len1, len2, delimiter;
		boolean dual=true;

		if(dual){
			len1=10;
			len2=10;
			delimiter='+';
			pairs.add(new Pair("ACGTACGTAC+CGAACGAACG","ACGTACGTAC+CGAACGAACG"));
			pairs.add(new Pair("ACGTACGTAC+NNNNNNNNNN","ACGTACGTAC+CGAACGAACG"));
			pairs.add(new Pair("GGGGGGGGGG+NNNNNNNNNN","GGGGGGGGGG+AGAACGAACG"));
			pairs.add(new Pair("GGGGGGGGGG+GNNNNNNNNN","GGGGGGGGGG+CGAACGAACG"));
			pairs.add(new Pair("GGGGGGGGGG+TCGTACGTTT","GGGGGGGGGG+CGAACGAACG"));
			pairs.add(new Pair("NNNNNNNNNN+AAAAAAAAAA","GGGGGTGGGG+AAAAAAAAAA"));
		}else{

			len1=10;
			len2=0;
			delimiter=0;
			pairs.add(new Pair("ACGTACGTAC","ACGTACGTAC"));
			pairs.add(new Pair("CCGTACGTAC","ACGTACGTAC"));
			pairs.add(new Pair("AGGGGGGGGG","GGGGGGGGGG"));
			pairs.add(new Pair("TGGGGGGGGG","GGGGGGGGGG"));
			pairs.add(new Pair("GGGGGGGGGG","GGGGGGGGGG"));
			pairs.add(new Pair("NNNNNNNNNN","GGGGGTGGGG"));
		}

		Collections.sort(pairs);
		System.err.println("Expected:\n"+pairs);
		System.err.println("\nWriting:");

		ByteArrayOutputStream baos=new ByteArrayOutputStream();
		DemuxServer.writeToStream(pairs, baos, coding, len1, len2);
		try{
			baos.close();
		}catch(IOException e){
// ByteArrayOutputStream.close() should not fail, but retain diagnostics.
			e.printStackTrace();
		}
		String result=baos.toString();
		System.err.println("Created:\n"+result);


		System.err.println("Decoding.");
		HashMap<String, String> map=null;
		if(coding==Sketch.RAW){
			map=toMapRaw(result, true);
		}else{
			map=toMapA48(result.getBytes(), true, coding, len1, len2, (byte)delimiter);
		}
		ArrayList<DemuxServer.Pair> pairs2=DemuxServer.toPairList(map);
		System.err.println("Decoded:\n"+pairs2);
	}

/** Creates a client using the shared default demux server URL. */
	public DemuxClient(){}

/**
	 * Creates a client with an optional server URL override.
	 * @param address_ Base demux server URL ending in '/', or null to use the
	 * shared default
	 */
	public DemuxClient(String address_){
		if(address_!=null){address=address_;}
	}

/**
	 * Sends demux data to the configured server and decodes its response.
	 *
	 * @param dd Demux data containing barcodes and configuration
	 * @param verbose Whether to print debugging information
	 * @return map of observed barcode keys to corrected barcode values
	 */
	HashMap<String, String> getMap(DemuxData dd, boolean verbose){
		String response=sendAndReceive(dd, verbose);
		if(dd.coding==Sketch.A48){
			return toMapA48(response.getBytes(), verbose, dd);
		}else{
			return toMapRaw(response, verbose);
		}
	}

/**
	 * Encodes one request, appends {@code demux} to the configured base URL,
	 * posts every chunk, and returns the raw response body.
	 * <p>
	 * If {@link #writeSent} is set, the exact request bytes are also attempted
	 * to that file for debugging; debug-write close errors are ignored.
	 * Non-2xx responses call {@link KillSwitch#kill()} unless
	 * {@link ServerTools#suppressErrors} is set; embedded legacy server errors
	 * kill regardless of that flag.
	 *
	 * @param dd Demux data to encode and send
	 * @param verbose Whether to print debugging information
	 * @return server response string
	 */
	private String sendAndReceive(DemuxData dd, boolean verbose){
		if(verbose && false){
			System.err.println("dd: type="+dd.type+", len1="+dd.length1+", len2="+dd.length2+
", delimiter="+dd.barcodeDelimiter+", hdistsum="+dd.hdistSum+
", expected="+dd.expectedList.size()+", counts="+dd.codeCounts.size());
		}

		ArrayList<byte[]> chunks=dd.encode();
		for(byte[] chunk : chunks){
			assert(chunk[chunk.length-1]=='\n') : "'"+new String(chunk, chunk.length-60, 60)+"'";
		}

		if(verbose){System.err.println("Sending "+chunks.size()+" chunks.");}

		if(writeSent!=null){
			ByteStreamWriter bsw=ByteStreamWriter.makeBSW(writeSent, true, false, false);
			for(byte[] chunk : chunks){
				bsw.print(chunk);
			}
//Best-effort debug copy; request success does not depend on this close.
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
		if(e){
			System.err.println("An error was encountered: \n"+response.s);
			KillSwitch.kill();
		}
		return response.s;
	}

/**
	 * Parses one raw {@code key<TAB>value} response pair per line.
	 *
	 * @param response Server response string in raw format
	 * @param verbose Whether to print debugging information
	 * @return map of parsed key-value pairs
	 */
	private static HashMap<String, String> toMapRaw(String response, boolean verbose){

		HashMap<String, String> map=new HashMap<String, String>();
		LineParserS2 responseParser=new LineParserS2('\n');
		LineParserS1 tp=new LineParserS1('\t');

		responseParser.set(response);
		while(responseParser.hasMore()){
			responseParser.advance();
			if(responseParser.currentFieldLength()>0){
//System.err.println(responseParser.a()+", "+responseParser.b());
				String line=responseParser.parseStringFromCurrentField();
				tp.set(line);
				String key=null, value=null;
				try {
					key=tp.parseString(0);
					value=tp.parseString(1);
				}catch(Throwable t){
					t.printStackTrace();
					System.err.println("'"+new String(line)+"'; "+tp.terms());
				}
				assert(!map.containsKey(key));
				map.put(key, value);
			}
		}

		return map;
	}

/**
	 * Parses A48-encoded response using DemuxData configuration.
	 *
	 * @param response Server response bytes in A48 format
	 * @param verbose Whether to print debugging information
	 * @param dd DemuxData containing encoding parameters
	 * @return map of decoded key-value pairs
	 */
	private static HashMap<String, String> toMapA48(byte[] response, boolean verbose, DemuxData dd){
		return toMapA48(response, verbose, dd.coding, dd.length1, dd.length2, (byte)dd.barcodeDelimiter);
	}

/**
	 * Parses one A48-encoded key/value response pair per line.
	 * <p>
	 * Server responses omit a field when it matches the previous line's field;
	 * this method restores omitted segments from {@code prev} before assembling
	 * the single-barcode or dual-barcode key/value strings.
	 *
	 * @param response Server response bytes in A48 format
	 * @param verbose Whether to print debugging information
	 * @param coding Encoding type retained with the response parameters; currently unused
	 * @param len1 Length of first barcode segment
	 * @param len2 Length of second barcode segment (0 for single mode)
	 * @param delimiter Character separating dual barcodes (0 for concatenation)
	 * @return map of decoded barcode pairs
	 */
	private static HashMap<String, String> toMapA48(byte[] response, boolean verbose, int coding, int len1, int len2, byte delimiter){
		HashMap<String, String> map=new HashMap<String, String>();
		LineParser2 responseParser=new LineParser2('\n');
		LineParser2 tp=new LineParser2('\t');
		final ByteBuilder bb=new ByteBuilder(32);

		byte[][] prev=new byte[4][];
		byte[][] current=new byte[4][];
		final int[] len=new int[] {len1, len2<1 ? len1 : len2, len1, len2};
		responseParser.set(response);
		while(responseParser.hasMore()){
			responseParser.advance();
			if(responseParser.currentFieldLength()>0){
//System.err.println(responseParser.a()+", "+responseParser.b());
				byte[] line=responseParser.parseByteArrayFromCurrentField();
				tp.set(line);
				int term=0;
/*
				 * DemuxServer.appendA48 emits an empty field to repeat the
				 * same segment from the previous pair. This decoder relies
				 * on the server always writing every delimiter: one tab in
				 * single-barcode mode and three in dual-barcode mode.
				 *
				 * A trailing empty field is load-bearing because Pair sorts
				 * by value first, so adjacent lines often share the value
				 * side. LineParser2 currently reports a final empty term by
				 * leaving b on the trailing delimiter after advance(). If
				 * that changed, the assertion below would trip on term=3.
				 */
				while(tp.hasMore()){
					tp.advance();
					if(tp.currentFieldLength()<1){
						current[term]=prev[term];
					}else{
						current[term]=decodeA48(line, tp.a(), tp.b(), bb, len[term]);
					}
					term++;
				}

				assert(term==2 || term==4) : term+": '"+new String(line)+"'";
				final String key, value;
				if(term==2){
					key=new String(current[0]);
					value=new String(current[1]);
				}else if(delimiter>0){
					key=bb.clear().append(current[0]).append(delimiter).append(current[1]).toString();
					value=bb.clear().append(current[2]).append(delimiter).append(current[3]).toString();
				}else{
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

/**
	 * Decodes a single A48-compressed barcode segment.
	 * Parses from {@code from} through the first byte below {@code '0'}, then
	 * converts the compressed numeric representation back to ACGTN bases.
	 *
	 * @param term byte array containing encoded data
	 * @param from start position in array, inclusive
	 * @param to unused retained field end from the tab parser
	 * @param bb ByteBuilder for constructing result
	 * @param len expected length of decoded barcode
	 * @return decoded barcode as byte array
	 */
	public static byte[] decodeA48(byte[] term, int from, int to, ByteBuilder bb, int len){
		bb.clear();
		long x=DemuxData.parseA48(term, from);
		DemuxData.decodeACGTN(x, bb, len);
		return bb.toBytes();
	}

/**
	 * Detects old server-side HTTP error text propagated through a response body.
	 * @param s Server response body, or null
	 * @return true if the body contains the historical HTTP error marker
	 */
	private static boolean checkForError(String s){
		if(s==null){return false;}
		return s.contains("java.io.IOException: Server returned HTTP response code:");
	}

/** True after a response body contains a server-side error marker. */
	public boolean errorState=false;
/** Base URL for demux HTTP requests; must end with '/' before sendAndReceive appends 'demux'. */
	public String address=Shared.demuxServer();
/** Optional debug output path for exact bytes sent to the server. */
	public static String writeSent=null;

}
