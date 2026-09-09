package ml;

import java.io.File;
import java.nio.charset.StandardCharsets;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/**
 * Gate runner: applies a textual BBNet to the data rows of a {@code #dims} TSV and
 * writes {@code row<TAB>out0[<TAB>out1...]} lines, one per input row.
 *
 * <p>Lives in package {@code ml} because the parser and the format contract are
 * owned here (CellNet's inference methods are public). It is a development and
 * measurement tool in the BBTools sense: it exists to check that an exported or
 * converted net reproduces a trainer's predictions row for row, not to score
 * production bins.</p>
 *
 * <p>Usage: {@code bbnetapply net=<file.bbnet> in=<vectors.tsv>
 * out=<predictions.tsv> rows=<n>}</p>
 *
 * <p>Failure contract (crash loud, never a silent partial): every malformed input
 * (missing or mismatched {@code #dims}, short row, non-finite input or output,
 * unknown or malformed flag, non-positive {@code rows=}, unreadable net/input file,
 * missing output directory) throws before or during the row loop;
 * an I/O error latched by the reader ({@link ByteFile#close()}) or the writer
 * ({@link ByteStreamWriter#poisonAndWait()}) throws AFTER both are closed, so an
 * exception thrown inside the loop is never masked by the close path; zero data
 * rows throws. Note that with BBTools 54b62510 a flush/close-time write failure is
 * still dropped UPSTREAM ({@code ReadWrite.finishWriting} discards {@code close()}'s
 * return; mag-qc {@code records/BBTOOLS_BUGS_FOUND_v1.md} row 13) - the one-line
 * upstream fix makes this runner reject it, as {@link BBNetApplyTest} shows when run
 * with {@code -Dmagqc.readwrite=patched}.</p>
 *
 * @author Sayu (2026-09-03), UMP45 (error propagation, 2026-09-08)
 */
public class BBNetApply {

	public static void main(String[] args){
		String netFile=null, inFile=null, outFile=null;
		int rows=Integer.MAX_VALUE;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Malformed argument (expected key=value): '"+arg+"'\n"+usage());}
			final String key=arg.substring(0, eq).toLowerCase();
			final String value=arg.substring(eq+1);
			if(key.equals("net")){netFile=value;}
			else if(key.equals("in")){inFile=value;}
			else if(key.equals("out")){outFile=value;}
			else if(key.equals("rows")){rows=Integer.parseInt(value);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg+"\n"+usage());}
		}
		if(netFile==null || inFile==null || outFile==null){
			throw new IllegalArgumentException(usage());
		}
		if(rows<=0){throw new IllegalArgumentException("rows must be positive, got "+rows+"\n"+usage());}
		//Preflight the two paths whose upstream failure is a hard JVM exit, not an exception
		//(ReadWrite.java:1345 exceptionKill on a missing input; the writer's open path likewise),
		//so a caller or a test sees an ordinary rejection instead of a dead JVM.
		if(!new File(netFile).canRead()){throw new IllegalArgumentException("net file not readable: "+netFile);}
		if(!new File(inFile).canRead()){throw new IllegalArgumentException("input TSV not readable: "+inFile);}
		final File outParent=new File(outFile).getAbsoluteFile().getParentFile();
		if(outParent==null || !outParent.isDirectory()){
			throw new IllegalArgumentException("output directory does not exist: "+outFile);
		}
		final CellNet net=CellNetParser.load(netFile, false);
		if(net==null){throw new IllegalArgumentException("Could not load net: "+netFile);}
		final ByteFile bf=ByteFile.makeByteFile(inFile, true);
		ByteStreamWriter bsw=null;
		int inputWidth=-1, headerOutputs=-1, written=0;
		boolean ioError=false;
		try{
			for(byte[] raw=bf.nextLine(); raw!=null; raw=bf.nextLine()){
				if(raw.length==0){continue;}
				final String line=new String(raw, StandardCharsets.UTF_8).trim();
				if(line.length()==0){continue;}
				if(line.startsWith("#dims")){
					final String[] fields=line.substring(5).trim().split("\\s+");
					if(fields.length<1){throw new IllegalArgumentException("Malformed #dims header: "+line);}
					inputWidth=Integer.parseInt(fields[0]);
					if(fields.length>1){headerOutputs=Integer.parseInt(fields[1]);}
					if(inputWidth!=net.numInputs()){
						throw new IllegalArgumentException("input width mismatch: TSV #dims="+inputWidth+
							" net="+net.numInputs());
					}
					if(headerOutputs>=0 && headerOutputs!=net.numOutputs()){
						throw new IllegalArgumentException("output width mismatch: TSV #dims="+headerOutputs+
							" net="+net.numOutputs());
					}
					continue;
				}
				if(line.charAt(0)=='#'){continue;}
				if(inputWidth<0){throw new IllegalArgumentException("Input TSV is missing #dims header");}
				if(written>=rows){break;}
				final String[] fields=line.split("\\t", -1);
				if(fields.length<inputWidth){
					throw new IllegalArgumentException("row "+written+" has "+fields.length+
						" columns; expected at least "+inputWidth);
				}
				final float[] input=new float[inputWidth];
				for(int i=0; i<input.length; i++){
					input[i]=Float.parseFloat(fields[i]);
					if(!Float.isFinite(input[i])){
						throw new IllegalArgumentException("row "+written+" input "+i+" is non-finite");
					}
				}
				net.applyInput(input);
				net.feedForward();
				if(bsw==null){
					//Opened at the first VALID data row, so a malformed header or first row never
					//leaves an output file behind.
					bsw=new ByteStreamWriter(outFile, true, false, true);
					bsw.start();
				}
				final ByteBuilder bb=new ByteBuilder(32+net.numOutputs()*16);
				bb.append(written);
				for(int k=0; k<net.numOutputs(); k++){
					final float output=net.getOutput(k);
					if(!Float.isFinite(output)){
						throw new IllegalStateException("row "+written+" output "+k+" is non-finite");
					}
					bb.tab().appendSlow(output);
				}
				bb.nl();
				bsw.print(bb);
				written++;
			}
		}finally{
			//Close both ends and LATCH their error flags. Nothing is thrown from here, so an
			//exception raised inside the loop propagates unmasked; the latched flags are
			//judged after the block. ByteFile1.close() returns its read-error latch
			//(ByteFile1.java:134-144); poisonAndWait() returns the writer's errorState
			//(ByteStreamWriter.java:281-285), which is set from ReadWrite.finishWriting
			//(:148) -- see the class javadoc for the upstream row-13 caveat.
			ioError=bf.close();
			if(bsw!=null){ioError=bsw.poisonAndWait()|ioError;}
		}
		if(inputWidth<0){throw new IllegalArgumentException("Input TSV is missing #dims header");}
		if(headerOutputs>=0 && headerOutputs!=net.numOutputs()){
			throw new IllegalArgumentException("output width mismatch: TSV #dims="+headerOutputs+
				" net="+net.numOutputs());
		}
		if(ioError){
			throw new RuntimeException("BBNetApply: I/O error reported by the input reader or the output writer;"+
				" in="+inFile+" out="+outFile+" -- the output must not be trusted");
		}
		if(written==0){
			throw new IllegalStateException("BBNetApply: zero data rows in "+inFile+"; nothing written to "+outFile);
		}
		System.err.println("BBNetApply: wrote "+written+" rows to "+outFile);
	}

	private static String usage(){
		return "Usage: bbnetapply net=<file.bbnet> in=<vectors.tsv> out=<predictions.tsv> rows=<n>";
	}
}
