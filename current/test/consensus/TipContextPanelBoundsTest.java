package test.consensus;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import consensus.TipContextPanel;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;

/** Tests native panel discovery bounds without changing evidence acceptance.
 * @author Fischl */
public final class TipContextPanelBoundsTest {

	/** Exercises short-overlap boundaries through the public CLI on both strands. */
	public static void main(final String[] args) throws Exception{
		check(args.length==0, "This deterministic fixture has no arguments.");
		for(boolean reverse : new boolean[]{false, true}){
			trial(33, -1, reverse, false, 33, 200, 1);
			trial(33, 16, reverse, false, null, 200, 0);
			trial(33, 16, reverse, false, 33, 200, 1);
			trial(33, 16, reverse, true, 33, 200, 2);
			trial(33, 15, reverse, false, 33, 200, 0);
			trial(59, 28, reverse, false, null, 200, 0);
			trial(59, 28, reverse, false, 33, 200, 1);
			trial(61, 29, reverse, false, null, 200, 1);
			trial(61, 29, reverse, false, 33, 60, 0);
		}
		boolean failed=false;
		try{TipContextPanel.main(new String[]{"mode=prepare", "minoverlap=32"});}
		catch(IllegalArgumentException expected){failed=expected.getMessage().contains("Overlap bounds");}
		check(failed, "A32bp minimum cannot contain both16bp anchors and a mismatch.");
		System.out.println("TIP_PANEL_BOUNDS_TEST_PASS trials=18 indexed_parity=true checks="+checks);
	}

	/** Writes an exact suffix/prefix overlap and verifies multiplicity and context emission. */
	private static void trial(final int overlap, final int mismatch, final boolean reverse,
			final boolean duplicate, final Integer minimum, final int maximum, final int expected) throws Exception{
		assert(overlap>mismatch && mismatch>=-1) : "A substitution lies inside the overlap; -1 selects an exact overlap.";
		final byte[] genome=new byte[600], alphabet={'A', 'C', 'G', 'T'};
		final Random random=new Random(252033);
		for(int i=0; i<genome.length; i++){genome[i]=alphabet[random.nextInt(4)];}
		byte[] target=Arrays.copyOfRange(genome, 0, 200);
		byte[] donor=Arrays.copyOfRange(genome, 200-overlap, 400);
		final int position=200-overlap+mismatch;
		if(mismatch>=0){target[position]=target[position]=='A' ? (byte)'C' : (byte)'A';}
		if(reverse){target=AminoAcid.reverseComplementBases(target); donor=AminoAcid.reverseComplementBases(donor);}
		final Path dir=Files.createTempDirectory("fischl-panel-bounds-");
		final String input=dir.resolve("assembly.fa").toString(), targets=dir.resolve("targets.tsv").toString();
		final String output=dir.resolve("panel").toString();
		final ByteBuilder fasta=new ByteBuilder(">target\n").append(target).nl().append(">donor\n").append(donor).nl();
		if(duplicate){fasta.append(">duplicate\n").append(donor).nl();}
		write(input, fasta);
		write(targets, new ByteBuilder("target\n"));
		final ArrayList<String> arguments=new ArrayList<String>();
		arguments.add("mode=prepare"); arguments.add("in="+input); arguments.add("targets="+targets);
		arguments.add("dir="+output); arguments.add("t=1"); arguments.add("maxoverlap="+maximum);
		if(minimum!=null){arguments.add("minoverlap="+minimum);}
		TipContextPanel.main(arguments.toArray(new String[0]));
		arguments.add("indexed=t");
		arguments.add("dir="+dir.resolve("indexed"));
		TipContextPanel.main(arguments.toArray(new String[0]));
		for(String name : new String[]{"targets.tsv", "candidates.tsv"}){
			check(Arrays.equals(Files.readAllBytes(dir.resolve("panel/"+name)), Files.readAllBytes(dir.resolve("indexed/"+name))),
					"Indexed discovery must preserve every candidate count, name, span and offset: "+name);
		}
		for(String side : new String[]{"left", "right"}){
			final String path="case001/"+side+".context.fa";
			if(Files.exists(dir.resolve("panel/"+path))){
				check(Arrays.equals(Files.readAllBytes(dir.resolve("panel/"+path)), Files.readAllBytes(dir.resolve("indexed/"+path))),
						"Indexed discovery must preserve the context sequence on "+side);
			}
		}
		final LineParser1 parser=new LineParser1('\t');
		int found=0;
		for(byte[] line : ByteFile.toLines(output+"/candidates.tsv")){
			parser.set(line);
			if(!parser.termEquals("case001", 0)){continue;}
			final boolean chosen=parser.termEquals(reverse ? "left" : "right", 1);
			check(parser.parseInt(2)==(chosen ? expected : 0), "Unexpected candidate count in "+output);
			if(chosen){
				found++;
				check(parser.parseInt(3)==(expected==1 ? overlap : -1), "Reported overlap differs from the fixture.");
			}
		}
		check(found==1, "Target must have exactly one recorded tested end.");
		check(Files.exists(dir.resolve("panel/case001/"+(reverse ? "left" : "right")+".context.fa"))==(expected==1),
				"Only a unique donor may produce a context.");
	}

	/** Uses the native writer and refuses silent partial fixture creation. */
	private static void write(final String path, final ByteBuilder bytes){
		assert(path!=null && bytes!=null) : "Fixture output needs a path and sequence bytes.";
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start(); writer.print(bytes);
		check(!writer.poisonAndWait(), "Failed to finish fixture "+path);
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}
	private static int checks;
}
