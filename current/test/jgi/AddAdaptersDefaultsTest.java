package test.jgi;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;

import jgi.AddAdapters;
import parse.Parse;
import stream.FASTQ;

/**
 * Constructor-only parsing-default and explicit-override controls for deprecated AddAdapters.
 * A fresh generated input and literal adapter let construction finish without calling write/read,
 * opening a main read stream or invoking the incompatible generic header-parsing path.
 * Baseline records the two omitted-flag defaults as true; fixed mode requires false.
 * Explicit aliases and last-option-wins expectations are identical in both phases.
 * @author Sucrose
 */
public final class AddAdaptersDefaultsTest{

	/** Requires {@code phase=baseline|fixed} and a new {@code outdir=...}; run in a fresh JVM. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown constructor-fixture argument: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null, "Set a new outdir for generated constructor inputs and observations.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		final String path=root.toString();
		require(path.indexOf('=')<0 && path.indexOf('\t')<0 && path.indexOf('\n')<0 && path.indexOf('\r')<0,
			"AddAdapters splits arguments at equals signs; fixture paths must also fit one table row.");
		Files.createDirectory(root);
		final Path input=root.resolve("ordinary.fq");
		write(input, INPUT);
		write(root.resolve("input-before.fq"), INPUT);
		write(root.resolve("cases.tsv"), "mode\tcase\tbefore\tintended\tphaseExpected\tactual\n");
		final boolean baseline=phase.equals("baseline");
		final boolean original=FASTQ.PARSE_CUSTOM;
		try{
			for(final String mode : new String[]{"write", "grade"}){
				for(int i=0; i<CASES.length; i++){
					final String[] row=CASES[i];
					final boolean intended=INTENDED[i];
					final boolean expected=i==0 ? baseline : intended;
					//Opposite inherited state proves the constructor/explicit options establish the observed value.
					FASTQ.PARSE_CUSTOM=!expected;
					final boolean before=FASTQ.PARSE_CUSTOM;
					final Path output=root.resolve("unused-"+mode+"-"+row[0]+".fq");
					final ArrayList<String> constructorArgs=new ArrayList<String>(Arrays.asList(
						"in="+input, "out="+output, "mode="+mode, "literal=AC", "adderrors=f", "arc=f",
						"int=f", "qin=33", "qout=33", "simd=f", "t=1", "ow=t", "app=f"));
					for(int j=1; j<row.length; j++){constructorArgs.add(row[j]);}
					Files.write(root.resolve(mode+"-"+row[0]+".argv"), constructorArgs,
						StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
					new AddAdapters(constructorArgs.toArray(new String[constructorArgs.size()]));
					final boolean actual=FASTQ.PARSE_CUSTOM;
					append(root.resolve("cases.tsv"), mode+"\t"+row[0]+"\t"+before+"\t"+
						intended+"\t"+expected+"\t"+actual+"\n");
					check(actual==expected, mode+"/"+row[0]+": expected phase parser flag="+expected+", observed="+actual);
					check(!Files.exists(output), mode+"/"+row[0]+": constructor must not create its unused main output.");
					check(matches(input, INPUT), mode+"/"+row[0]+": constructor changed the generated input bytes.");
					if(actual!=intended){defaultMismatches++;}
					cases++;
				}
			}
		}finally{
			FASTQ.PARSE_CUSTOM=original;
		}
		Files.write(root.resolve("input-after.fq"), Files.readAllBytes(input), StandardOpenOption.CREATE_NEW);
		check(FASTQ.PARSE_CUSTOM==original, "Fixture must restore its inherited parser flag after constructor observations.");
		check(cases==12, "Two modes and six option sequences require exactly12 constructor cases: "+cases);
		check(defaultMismatches==(baseline ? 2 : 0),
			"Only the omitted-flag default differs from intended behavior in each mode before repair: "+defaultMismatches);
		if(failures>0){throw new AssertionError("AddAdapters constructor fixture failed "+failures+" of "+checks+" checks; inspect "+root);}
		System.out.println("ADD_ADAPTERS_DEFAULTS_PASS phase="+phase+" cases="+cases+
			" defaultMismatches="+defaultMismatches+" checks="+checks);
	}

	/** Exact byte comparison remains independent of FASTQ parsing and its global flags. */
	private static boolean matches(final Path path, final String expected) throws Exception{
		return Files.isRegularFile(path) && Arrays.equals(Files.readAllBytes(path), expected.getBytes(StandardCharsets.US_ASCII));
	}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}
	private static void append(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}
	}

	/** Label followed by literal options; aliases appear after earlier options to check their normal precedence. */
	private static final String[][] CASES={
		{"omitted"},
		{"explicit-f", "parsecustom=f"},
		{"explicit-t", "parsecustom=t"},
		{"alias-t", "fastqparsecustom=t"},
		{"t-then-alias-f", "parsecustom=t", "fastqparsecustom=f"},
		{"f-then-alias-t", "parsecustom=f", "fastqparsecustom=t"}
	};
	private static final boolean[] INTENDED={false, false, true, true, false, true};
	private static final String INPUT="@ordinary\nG\n+\nI\n";
	private static int checks, failures, cases, defaultMismatches;
}
