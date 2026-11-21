package test.stream;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import shared.Shared;
import stream.Read;
import stream.ReadStreamSamWriter;
import stream.ReadStreamWriter;
import stream.SamLine;

/** Exercises normal SAM adapter completion or deliberately invalid forwarding state.
 * The failure invocation must halt with the adapter's diagnostic and original stack;
 * a successful return or merely ending the adapter thread is not acceptance.
 * Use a fresh JVM with assertions enabled and a new output file for each invocation.
 * @author Shinobu
 */
public final class SamAdapterFailureTest{

	/** Requires out=file.sam and failure=t|f; invalid forwarding is never a baseline liveness run. */
	public static void main(final String[] args) throws Exception{
		String output=null;
		Boolean failure=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("out")){output=pair[1];}
			else if(pair[0].equals("failure")){failure=Parse.parseBoolean(pair[1]);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(output==null || failure==null){throw new IllegalArgumentException("Set out and failure.");}
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("This control requires assertions enabled.");}
		if(Data.scaffoldLocs!=null){throw new IllegalStateException("This control requires a fresh JVM without a reference.");}
		final Path path=Paths.get(output);
		if(Files.exists(path)){throw new IllegalArgumentException("Output already exists: "+path);}
		Shared.setThreads(1);
		ReadWrite.USE_READ_STREAM_SAM_WRITER=true;
		ReadStreamWriter.USE_ATTACHED_SAMLINE=!failure;
		SamLine.SET_FROM_OK=!failure;
		final String record="r0\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII";
		final SamLine sam=new SamLine(new LineParser1('\t').set(record.getBytes(StandardCharsets.US_ASCII)));
		final Read read=new Read("ACGT".getBytes(StandardCharsets.US_ASCII), new byte[]{40, 40, 40, 40}, "r0", 0);
		read.samline=sam;
		final FileFormat format=FileFormat.testOutput(output, FileFormat.SAM, null, false, false, false, false);
		final ReadStreamSamWriter writer=new ReadStreamSamWriter(format, 1, "@HD\tVN:1.6\tSO:unknown", false);
		writer.start();
		//A capacity-one outer queue also permits a producer to block behind the failed consumer.
		//The exact race is not asserted; fatal termination must work regardless of that schedule.
		for(int i=0; i<(failure ? 32 : 1); i++){
			final ArrayList<Read> batch=new ArrayList<Read>(1);
			batch.add(read);
			writer.addList(batch);
		}
		writer.poison();
		writer.join();
		if(failure){throw new AssertionError("The invalid forwarding state must terminate the JVM from the adapter.");}
		if(writer.errorState() || !writer.finishedSuccessfully()){
			throw new AssertionError("Normal adapter completion reported failure.");
		}
		int records=0;
		for(final String line : Files.readAllLines(path, StandardCharsets.US_ASCII)){
			if(!line.startsWith("@") && !line.isEmpty()){
				if(!line.equals(record)){throw new AssertionError("Unexpected SAM record: "+line);}
				records++;
			}
		}
		if(records!=1){throw new AssertionError("Expected one complete SAM record, observed "+records);}
		System.out.println("SAM_ADAPTER_NORMAL_PASS records=1");
	}
}
