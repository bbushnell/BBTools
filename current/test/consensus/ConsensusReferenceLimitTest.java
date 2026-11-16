package test.consensus;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import fileIO.ReadWrite;
import parse.Parse;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;

/** Checks that an alignment quota does not truncate any tool's reference input.
 * Invokes the actual private reference-reader factory, without starting mapping
 * workers that could wait on an absent reference. Runs in its own test JVM.
 * @author Shinobu */
public final class ConsensusReferenceLimitTest{

	/** Requires phase=baseline|fixed and a new outdir with an existing parent. */
	public static void main(final String[] args) throws Exception{
		String phase=null, outdir=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){outdir=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		if(!"baseline".equals(phase) && !"fixed".equals(phase)){
			throw new IllegalArgumentException("Set phase=baseline or phase=fixed.");
		}
		if(outdir==null){throw new IllegalArgumentException("Set outdir to a new directory.");}
		final Path root=Files.createDirectory(Paths.get(outdir)).toAbsolutePath();
		final Path ref=root.resolve("references.fa"), input=root.resolve("empty.sam");
		Files.write(ref, ">first\nACGTACGT\n>second\nTTTTCCCC\n".getBytes(StandardCharsets.US_ASCII));
		Files.write(input, ("@HD\tVN:1.6\n@SQ\tSN:first\tLN:8\n"+
			"@SQ\tSN:second\tLN:8\n").getBytes(StandardCharsets.US_ASCII));
		final StringBuilder observations=new StringBuilder("tool\talignment_limit\treferences\n");
		for(final String name : new String[]{"ConsensusMaker", "FixScaffoldGaps", "Lilypad"}){
			for(final long limit : new long[]{1, -1}){
				final List<String> expected="baseline".equals(phase) && limit==1 ?
					Arrays.asList("first") : Arrays.asList("first", "second");
				final Class<?> type=Class.forName("consensus."+name);
				final String[] options={"in="+input, "ref="+ref, "out="+root.resolve(name+limit+".fa"),
					"reads="+limit, "t=1", "simd=f", "ow=f"};
				final Object tool=type.getConstructor(String[].class).newInstance((Object)options);
				final Method factory=type.getDeclaredMethod("makeRefCris");
				factory.setAccessible(true);
				final ConcurrentReadInputStream stream;
				try{stream=(ConcurrentReadInputStream)factory.invoke(tool);}
				catch(InvocationTargetException e){throw new AssertionError("Reference factory failed for "+name, e.getCause());}
				final List<String> names=new ArrayList<String>();
				boolean failed;
				try{
					for(ListNum<Read> batch=stream.nextList(); batch!=null; batch=stream.nextList()){
						final boolean terminal=batch.isEmpty();
						for(final Read read : batch){names.add(read.id);}
						stream.returnList(batch);
						if(terminal){break;}
					}
				}finally{failed=ReadWrite.closeStream(stream);}
				observations.append(name).append('\t').append(limit).append('\t').append(names).append('\n');
				Files.write(root.resolve("observations.tsv"), observations.toString().getBytes(StandardCharsets.US_ASCII));
				if(failed || !names.equals(expected)){
					throw new AssertionError(name+" reads="+limit+" reference names="+names+", expected="+expected+", inputError="+failed);
				}
			}
		}
		System.out.println("CONSENSUS_REFERENCE_LIMIT_PASS phase="+phase+" tools=3 cases=6");
	}
}
