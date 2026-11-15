package test.hiseq;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;

import hiseq.MicroTile;
import stream.Read;

/** Generates finite primary/auxiliary quality fixtures and checks the model API.
 * CLI failures run in separate JVMs through the committed validation controller.
 * @author Brian Bushnell */
public class HiseqQualityBoundaryTest{

	/** Requires phase=baseline|fixed and a new outdir; never launches child programs. */
	public static void main(String[] args) throws Exception{
		String phase=null, directory=null;
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			require(pair.length==2, "Use name=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "phase=baseline|fixed required");
		require(directory!=null && new File(directory).mkdir(), "Use a new fixture directory");
		final boolean fixed=phase.equals("fixed");
		final String name="INST:1:FC:1:1101:100:100";
		final String bases="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
		final byte[] quality=new byte[bases.length()];
		Arrays.fill(quality, (byte)'?');
		final String qualities=new String(quality, StandardCharsets.US_ASCII);
		write(directory, "valid.fq", "@"+name+" 1:N:0:ACGT\n"+bases+"\n+\n"+qualities+"\n");
		write(directory, "missing.fa", ">"+name+"\n"+bases+"\n");
		write(directory, "missing.sam", name+"\t4\t*\t0\t0\t*\t*\t0\t0\t"+bases+"\t*\n");
		write(directory, "missing-mate.sam", name+"\t77\t*\t0\t0\t*\t*\t0\t0\t"+bases+"\t"+qualities+"\n"
			+name+"\t141\t*\t0\t0\t*\t*\t0\t0\t"+bases+"\t*\n");
		write(directory, "auxiliary.sam", "@SQ\tSN:chr1\tLN:1000\n"+name+"\t0\tchr1\t1\t60\t"
			+bases.length()+"=\t*\t0\t0\t"+bases+"\t*\n");
		write(directory, "empty.fq", "");
		System.out.println("case\texception\treads\tbases\tresult");
		probe(false, fixed);
		probe(true, fixed);
	}

	/** Checks rejection before model mutation, alongside a real-quality control. */
	private static void probe(boolean quick, boolean fixed){
		final MicroTile missing=new MicroTile();
		final Read read=new Read(new byte[] {'A', 'C', 'G', 'T'}, null, "missing-quality", 0);
		RuntimeException caught=null;
		try{if(quick){missing.addQuick(read);}else{missing.add(read);}}
		catch(RuntimeException e){caught=e;}
		require(caught!=null, "Missing qualities must not be silently accepted");
		require(fixed ? caught instanceof IllegalArgumentException : caught instanceof NullPointerException,
			"Unexpected missing-quality exception: "+caught);
		require(missing.readCount==(fixed ? 0 : 1) && missing.baseCount==(fixed ? 0 : 4),
			"Rejected reads must not mutate model counters: "+missing.readCount+", "+missing.baseCount);
		if(fixed){require(caught.getMessage()!=null && caught.getMessage().toLowerCase().contains("qualit"),
			"Diagnostic must explain missing qualities: "+caught);}
		System.out.println((quick ? "addQuick" : "add")+"\t"+caught.getClass().getSimpleName()+"\t"
			+missing.readCount+"\t"+missing.baseCount+"\tPASS");
		final MicroTile control=new MicroTile();
		final Read good=new Read(new byte[] {'A', 'C', 'G', 'T'}, new byte[] {30, 30, 30, 30}, "quality-control", 1);
		if(quick){control.addQuick(good);}else{control.add(good);}
		require(control.readCount==1 && control.baseCount==4, "Valid qualities must retain the read");
		control.add((Read)null);
		control.addQuick(null);
		control.add(new Read(new byte[0], null, "empty-control", 2));
		require(control.readCount==1 && control.baseCount==4, "Null/empty reads must remain no-ops");
		System.out.println((quick ? "quick-control" : "control")+"\t-\t1\t4\tPASS");
	}

	/** Writes only small literal fixture inputs in the requested new directory. */
	private static void write(String directory, String name, String text) throws Exception{
		Files.write(new File(directory, name).toPath(), text.getBytes(StandardCharsets.US_ASCII));
	}

	/** Keeps fixture checks enabled independently of JVM assertion settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
