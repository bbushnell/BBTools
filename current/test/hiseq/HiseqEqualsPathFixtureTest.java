package test.hiseq;

import java.io.File;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Locale;

import hiseq.BGI2Illumina;
import hiseq.PlotFlowCell;
import hiseq.PlotReadPosition;

/** Finite hiseq equals-in-file-argument parser checks.
 * @author Brian Bushnell, Jean */
public class HiseqEqualsPathFixtureTest{

	/** Accepts phase=baseline|fixed, scenario=bgi|plotflowcell|plotreadposition, outdir=fresh. */
	public static void main(String[] args) throws Exception{
		String phase=null, scenario=null, directory=null;
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			require(pair.length==2, "Use name=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("scenario")){scenario=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "phase=baseline|fixed required");
		require(directory!=null, "outdir=fresh is required");
		require(scenario!=null, "scenario is required");
		Locale.setDefault(Locale.ROOT);
		final File out=fresh(directory);
		if(scenario.equals("bgi")){
			bgi(out, phase.equals("fixed"));
		}else if(scenario.equals("plotflowcell")){
			plotFlowCell(out, phase.equals("fixed"));
		}else if(scenario.equals("plotreadposition")){
			plotReadPosition(out, phase.equals("fixed"));
		}else{throw new IllegalArgumentException("Unknown scenario: "+scenario);}
	}

	private static void bgi(File out, boolean fixed) throws Exception{
		final Paths paths=paths(out);
		final BGI2Illumina tool=new BGI2Illumina(new String[] {
				"in="+paths.literal,
				"out="+new File(out, "bgi.fq").getPath(),
				"parseextra=f"});
		check("BGI2Illumina.in1", stringField(tool, "in1"), paths, fixed);
	}

	private static void plotFlowCell(File out, boolean fixed) throws Exception{
		final Paths paths=paths(out);
		final PlotFlowCell tool=new PlotFlowCell(new String[] {
				"in="+paths.literal,
				"dump="+new File(out, "flowcell.tsv").getPath(),
				"loadkmers=f"});
		check("PlotFlowCell.in1", stringField(tool, "in1"), paths, fixed);
	}

	private static void plotReadPosition(File out, boolean fixed) throws Exception{
		final Paths paths=paths(out);
		final Paths barcodes=barcodes(out);
		final PlotReadPosition tool=new PlotReadPosition(new String[] {
				"in="+paths.literal,
				"expected="+barcodes.literal,
				"out="+new File(out, "positions.tsv").getPath()});
		check("PlotReadPosition.in1", stringField(tool, "in1"), paths, fixed);
		check("PlotReadPosition.expectedPath", stringField(tool, "expectedPath"), barcodes, fixed);
	}

	private static void check(String label, String observed, Paths paths, boolean fixed){
		final String expected=fixed ? paths.literal : paths.prefix;
		require(expected.equals(observed),
				label+" expected "+expected+" but observed "+observed);
		System.out.println(label+"\tPASS\t"+observed);
	}

	private static Paths paths(File out) throws Exception{
		final File prefix=new File(out, "input.fq");
		final File literal=new File(out, "input.fq=literal.fq");
		write(prefix, "@prefix:1:FC:1:1101:1:1 1:N:0:AAAA\nACGT\n+\nIIII\n");
		write(literal, "@literal:1:FC:1:1101:1:1 1:N:0:CCCC\nTGCA\n+\nIIII\n");
		return new Paths(prefix, literal);
	}

	private static Paths barcodes(File out) throws Exception{
		final File prefix=new File(out, "barcodes.txt");
		final File literal=new File(out, "barcodes.txt=literal.txt");
		write(prefix, "AAAA\n");
		write(literal, "CCCC\n");
		return new Paths(prefix, literal);
	}

	private static void write(File file, String text) throws Exception{
		Files.write(file.toPath(), text.getBytes(StandardCharsets.US_ASCII));
	}

	private static String stringField(Object object, String field) throws Exception{
		Field f=object.getClass().getDeclaredField(field);
		f.setAccessible(true);
		return (String)f.get(object);
	}

	private static File fresh(String directory){
		final File out=new File(directory);
		require(out.isAbsolute() && out.mkdir(), "Cannot create fresh outdir: "+out);
		return out;
	}

	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static class Paths{
		Paths(File prefix_, File literal_){
			prefix=prefix_.getPath();
			literal=literal_.getPath();
		}
		final String prefix;
		final String literal;
	}
}
