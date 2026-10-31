package test.json;

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;

import json.JsonObject;
import json.JsonParser;
import parse.Parse;
import structures.ByteBuilder;

/** Measures ordinary ASCII JSON parsing or serialization with retained results.
 * @author Shinobu */
public final class JsonBenchmark{

	/** Runs one workload; use external wall/CPU timing for comparisons.
	 * @param args mode=parse|write and operations=positive count
	 * @throws Exception If the SHA-256 provider is unavailable */
	public static void main(final String[] args) throws Exception{
		String mode="parse";
		long operations=100000;
		for(String arg : args){
			final int split=arg.indexOf('=');
			if(split<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=arg.substring(0, split), value=arg.substring(split+1);
			if(key.equalsIgnoreCase("mode")){mode=value;}
			else if(key.equalsIgnoreCase("operations")){operations=Parse.parseKMG(value);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(operations<=0 || !(mode.equals("parse") || mode.equals("write"))){
			throw new IllegalArgumentException("Expected mode=parse|write and operations>0");
		}
		final byte[] input=TEXT.getBytes(StandardCharsets.US_ASCII);
		final JsonObject object=new JsonParser(input).parseJsonObject();
		if(object.getLong("count")!=314159L || object.getArray("hits").length!=2){
			throw new AssertionError("Benchmark input must retain count314159 and two hits");
		}
		final byte[] output=object.toText().toBytes();
		long checksum=0;
		Object last=null;
		final long start=System.nanoTime();
		if(mode.equals("parse")){
			for(long i=0; i<operations; i++){
				final JsonObject parsed=new JsonParser(input).parseJsonObject();
				checksum+=parsed.getLong("count")+parsed.getArray("hits").length;
				last=parsed;
			}
		}else{
			for(long i=0; i<operations; i++){
				final ByteBuilder bytes=object.toText();
				checksum+=bytes.length+(bytes.array[(int)(i%bytes.length)]&255);
				last=bytes;
			}
		}
		final long elapsed=System.nanoTime()-start;
		retained=last;
		System.out.println("JSON_BENCHMARK mode="+mode+" operations="+operations+
				" checksum="+checksum+" compute_ns="+elapsed+
				" input_sha80="+sha80(input)+" output_sha80="+sha80(output));
	}

	/** Produces the final 80 bits of SHA-256 for fixed workload identity.
	 * @param bytes Input bytes
	 * @return Twenty lowercase hexadecimal digits
	 * @throws Exception If the SHA-256 provider is unavailable */
	private static String sha80(final byte[] bytes) throws Exception{
		final byte[] digest=MessageDigest.getInstance("SHA-256").digest(bytes);
		final char[] out=new char[20], hex="0123456789abcdef".toCharArray();
		for(int i=0; i<10; i++){
			final int value=digest[digest.length-10+i]&255;
			out[i*2]=hex[value>>>4]; out[i*2+1]=hex[value&15];
		}
		return new String(out);
	}

	/** Prevents construction of this standalone benchmark. */
	private JsonBenchmark(){}
	/** Keeps the final generated object observable after the timed loop. */
	private static volatile Object retained;
	/** A fixed reply containing strings, numbers, flags, arrays and nested hits. */
	private static final String TEXT="{\"name\":\"sample\\tA\",\"count\":314159,\"score\":0.875,"
			+"\"flags\":[true,false,null],\"hits\":[{\"name\":\"alpha\",\"id\":7},"
			+"{\"name\":\"beta\",\"id\":9}]}";
}
