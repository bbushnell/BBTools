package test.barcode;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Locale;

import barcode.Barcode;
import barcode.PCRMatrix;
import barcode.PCRMatrixHDist;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Qualifies assignment-map write failure and checks nested instance-lock release.
 * The failing case requires the existing Unix character device /dev/full.
 * @author Shinobu
 */
public final class MatrixMapWriterFailureTest{

	/** Requires mode=baseline|fixed, api=static|instance, case=ok|full, and a new dir=path. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, MatrixMapWriterFailureTest.class, false);
		try{
			String mode=null, api=null, testCase=null, directory=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("api")){api=value;}
				else if(key.equals("case")){testCase=value;}
				else if(key.equals("dir")){directory=value;}
				else if(!Parser.parseStatic(arg, key, value)){
					throw new IllegalArgumentException("Unknown argument: "+arg);
				}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or mode=fixed.");
			check("static".equals(api) || "instance".equals(api), "Set api=static or api=instance.");
			check("ok".equals(testCase) || "full".equals(testCase), "Set case=ok or case=full.");
			check(directory!=null, "Set a new scratch directory.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Map-writer checks require -ea.");
			final Path dir=Files.createDirectory(Paths.get(directory));
			final boolean failure="full".equals(testCase);
			final Path target=failure ? Paths.get("/dev/full") : dir.resolve("map.tsv");
			if(failure){
				final int unixMode=((Number)Files.getAttribute(target, "unix:mode")).intValue();
				check((unixMode&0170000)==0020000 && Files.isWritable(target),
						"The failing case requires an existing writable character device /dev/full.");
			}

			final HashMap<String, String> assignments=new HashMap<String, String>();
			assignments.put("CCCC", "CCCC");
			assignments.put("AAAT", "AAAA");
			final HashMap<String, Barcode> counts=new HashMap<String, Barcode>();
			counts.put("AAAT", new Barcode("AAAT", 5));
			final PCRMatrixHDist matrix="instance".equals(api) ? new PCRMatrixHDist(4, 0, 0, false) : null;
			final ByteArrayOutputStream captured=new ByteArrayOutputStream();
			final PrintStream originalErr=System.err;
			Throwable thrown=null;
			try(final PrintStream capture=new PrintStream(captured, true, "UTF-8")){
				System.setErr(capture);
				try{
					if(matrix==null){
						PCRMatrix.printAssignmentMapStatic(assignments, target.toString(), counts, true, false);
					}else{
						//The collection overload delegates through the second, independently locking overload.
						matrix.printAssignmentMap(assignments, target.toString(), counts.values(), true, false);
					}
				}catch(final Throwable t){thrown=t;}
				finally{System.setErr(originalErr);}
			}
			final String diagnostics=captured.toString("UTF-8");
			originalErr.print(diagnostics);
			final boolean writerFailure=diagnostics.contains("java.io.IOException") &&
					diagnostics.contains("fileIO.ByteStreamWriter.run");
			if(failure){
				check(writerFailure, "Expected an actual asynchronous writer IOException; stderr:\n"+diagnostics);
				if("baseline".equals(mode)){
					check(thrown==null, "Baseline must return despite writer failure; thrown="+thrown);
				}else{
					check(thrown instanceof IllegalStateException &&
							String.valueOf(thrown.getMessage()).contains("assignment map") &&
							String.valueOf(thrown.getMessage()).contains("/dev/full"),
							"Fixed map writer must identify the failed assignment-map destination; thrown="+thrown);
				}
			}else{
				check(thrown==null && !writerFailure, "Valid map output failed: "+thrown+"\n"+diagnostics);
				final String actual=new String(Files.readAllBytes(target), StandardCharsets.US_ASCII);
				check("AAAT\tAAAA\t5\t1\nCCCC\tCCCC\t0\t0\n".equals(actual),
						"Sorted map rows/counts/distances changed: "+actual);
			}
			if(matrix!=null){
				final boolean acquired=matrix.rwlock().writeLock().tryLock();
				if(acquired){matrix.rwlock().writeLock().unlock();}
				check(acquired, "Printing must release both nested read locks before returning or throwing.");
			}
			System.out.println("MATRIX_MAP_WRITER mode="+mode+" api="+api+" case="+testCase+
					" writerFailure="+writerFailure+" propagated="+(thrown!=null)+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Always enforces the exact output/failure contracts, independently of downstream assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private MatrixMapWriterFailureTest(){}
}
