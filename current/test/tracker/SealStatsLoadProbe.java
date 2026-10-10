package test.tracker;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPOutputStream;

import fileIO.ByteFile1;
import fileIO.ReadWrite;
import shared.Shared;
import tracker.SealStats;

/** Finite actual-loader JT013 observations; no resource-leak measurement.
 * Each invocation uses one short input, built with Java only.
 * @author Brian Bushnell, Jean */
public final class SealStatsLoadProbe {

	public static void main(String[] args) throws IOException{
		assert(args.length==2 && args[0].startsWith("case=") && args[1].startsWith("out=")) :
			"JT013 takes one literal case and a fresh output directory.";
		final String name=args[0].substring(5);
		final Path out=Paths.get(args[1].substring(4));
		final boolean malformed=name.equals("malformed");
		final boolean corrupt=name.equals("trailer");
		final boolean gzip=corrupt || name.equals("valid_gzip");
		if(!malformed && !gzip && !name.equals("valid")){
			throw new IllegalArgumentException("Unknown JT013 fixture: "+name);
		}
		assert(!Files.exists(out)) : "Preserve earlier JT013 evidence: "+out;
		Files.createDirectories(out);
		final String content="#File\tsample.fq\n#Total\t5\t50\n#Matched\t5\t100%\n"+
			"alpha\t"+(malformed ? "not_a_count" : "5")+"\t100%\t50\t100%\n";
		byte[] bytes=content.getBytes(StandardCharsets.US_ASCII);
		if(gzip){bytes=gzip(bytes);}
		if(corrupt){
			assert(bytes.length>18) : "A gzip header and eight-byte trailer must surround the literal data.";
			bytes[bytes.length-8]^=1;//Corrupt only the first CRC32 trailer byte.
		}
		final Path input=out.resolve(gzip ? "input.tsv.gz" : "input.tsv");
		Files.write(input, bytes);
		Shared.setThreads(1);//ReadWrite.getGZipInputStream takes Java GZIPInputStream at one thread.
		ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_UNBGZIP=false;
		ReadWrite.USE_GUNZIP=false;
		ByteFile1.verbose=true;//Observe close/error diagnostics, not an OS resource-leak metric.
		System.err.println("JT013 case="+name+" threads="+Shared.threads()+" SIMD="+Shared.SIMD+" java_gzip=true");
		final String expected=(malformed || corrupt ? "THREW" : "RETURNED:5,50,1");
		String observed;
		try{
			final SealStats stats=new SealStats(input.toString());
			observed="RETURNED:"+stats.totalReads+","+stats.matchedBases+","+stats.map.size();
		}catch(RuntimeException e){
			observed="THREW:"+e.getClass().getName();
			e.printStackTrace(System.err);
		}catch(AssertionError e){
			observed="THREW:"+e.getClass().getName();
			e.printStackTrace(System.err);
		}
		final boolean agrees=expected.equals("THREW") ? observed.startsWith("THREW:") : expected.equals(observed);
		System.out.println(name+"\t"+expected+"\t"+observed+"\t"+agrees);
	}

	private static byte[] gzip(byte[] plain) throws IOException{
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(GZIPOutputStream gzip=new GZIPOutputStream(bytes)){gzip.write(plain);}
		return bytes.toByteArray();
	}
}
