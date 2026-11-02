package test.repeat;

import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import repeat.RepeatFinder;

/** Exercises real file parsing and the RepeatFinder command-line processing path.
 * The bounded contract is explicit rejection of paired input, not support for
 * scanning mates. Run each case in a fresh JVM. The unpaired control reads both
 * records separately and must report the second record's eight-base entropy region.
 * Linux descriptor and reader-stack checks verify closure after return/rejection.
 * @author Sucrose
 */
public final class RepeatPairedInputTest{

	/** Requires case=in2|forced|auto|unpaired and work=fresh-directory. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path work=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){work=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		require("in2".equals(mode) || "forced".equals(mode) || "auto".equals(mode) || "unpaired".equals(mode),
				"Choose case=in2, forced, auto or unpaired");
		require(work!=null, "work=fresh-directory is required");
		Files.createDirectory(work);
		final Path first=work.resolve("first.fq"), second=work.resolve("second.fq"), both=work.resolve("both.fq");
		final Path output=work.resolve("report.tsv");
		Files.write(first, FIRST.getBytes(StandardCharsets.US_ASCII));
		Files.write(second, SECOND.getBytes(StandardCharsets.US_ASCII));
		Files.write(both, (FIRST+SECOND).getBytes(StandardCharsets.US_ASCII));
		final ArrayList<String> command=new ArrayList<String>();
		command.add("in="+(mode.equals("in2") ? first : both));
		if(mode.equals("in2")){command.add("in2="+second);}
		command.add("interleaved="+(mode.equals("forced") ? "t" : mode.equals("auto") ? "auto" : "f"));
		command.add("out="+output);
		command.add("depth=f");
		command.add("entropy=t");
		command.add("ek=1");
		command.add("ew=4");
		command.add("minrepeat=4");
		command.add("k=3");
		command.add("mask=f");
		command.add("preview=27");
		command.add("qin=33");
		command.add("t=1");
		IllegalArgumentException rejection=null;
		try{RepeatFinder.main(command.toArray(new String[0]));}
		catch(IllegalArgumentException failure){rejection=failure;}
		final int readers=waitForReaders();
		final int descriptors=openInputs(first, second, both);
		System.out.println("case="+mode+" rejection="+rejection+" liveReaders="+readers+" openInputs="+descriptors);
		require(readers==0 && descriptors==0, "Returning or rejecting must not leave an owned reader open");
		if(mode.equals("unpaired")){
			require(rejection==null, "Explicit unpaired input was rejected: "+rejection);
			final List<String> rows=Files.readAllLines(output, StandardCharsets.US_ASCII);
			require(rows.size()==2, "Expected one header and one region from the second unpaired sequence");
			final String[] row=rows.get(1).split("\t", -1);
			require(row.length==14 && row[1].equals("8") && row[5].equals("1")
					&& row[6].equals("0") && row[7].equals("7") && row[12].equals("sample/2")
					&& row[13].equals("AAAAAAAA"), "Wrong unpaired region: "+rows.get(1));
		}else{
			require(rejection!=null && rejection.getMessage()!=null
					&& rejection.getMessage().toLowerCase(java.util.Locale.ROOT).contains("paired"),
					"Unsupported paired input must receive a clear IllegalArgumentException");
			require(!Files.exists(output), "Rejected input must not create a repeat report");
		}
		System.out.println("RepeatPairedInputTest "+mode+" PASS");
	}

	/** Allows brief thread termination bookkeeping; no process is killed or interrupted. */
	private static int waitForReaders() throws InterruptedException{
		final long end=System.nanoTime()+2000000000L;
		int live;
		while((live=liveReaders())>0 && System.nanoTime()<end){Thread.sleep(10);}
		return live;
	}

	/** Finds live reader coordinator/backend frames, including anonymous Thread wrappers. */
	private static int liveReaders(){
		int count=0;
		for(Map.Entry<Thread, StackTraceElement[]> entry : Thread.getAllStackTraces().entrySet()){
			if(entry.getKey()==Thread.currentThread() || !entry.getKey().isAlive()){continue;}
			for(StackTraceElement frame : entry.getValue()){
				final String name=frame.getClassName();
				if(name.startsWith("stream.ConcurrentGenericReadInputStream") || name.startsWith("fileIO.ByteFile2")){
					count++;
					break;
				}
			}
		}
		return count;
	}

	/** Counts descriptors still referring to the three owned fixture files on Linux. */
	private static int openInputs(Path... inputs) throws Exception{
		final Path descriptors=Paths.get("/proc/self/fd");
		require(Files.isDirectory(descriptors), "Linux /proc/self/fd is required for the closure check");
		int count=0;
		try(DirectoryStream<Path> entries=Files.newDirectoryStream(descriptors)){
			for(Path entry : entries){
				final Path target;
				try{target=Files.readSymbolicLink(entry);}
				catch(java.nio.file.NoSuchFileException closed){continue;}
				for(Path input : inputs){
					if(target.equals(input.toAbsolutePath().normalize())){count++;break;}
				}
			}
		}
		return count;
	}

	/** Explicit checks remain enabled with JVM assertions disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Only the second record contains a low-entropy region; paired names permit auto-detection. */
	private static final String FIRST="@sample/1\nACGTACGT\n+\nIIIIIIII\n";
	private static final String SECOND="@sample/2\nAAAAAAAA\n+\nIIIIIIII\n";
}
