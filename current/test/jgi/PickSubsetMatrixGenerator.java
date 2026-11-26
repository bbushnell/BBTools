package test.jgi;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Locale;
import java.util.Random;

import fileIO.FileFormat;
import parse.Parse;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Generates a deterministic complete numeric matrix for PickSubset timing checks.
 * This synthetic directed matrix is not a sequence-derived similarity model.
 * Includes one row per ordered pair, including self rows, in the twelve-column layout.
 * Uses bounded row storage and cached node names; output must be a new raw file.
 * Resolves the BBTools filename first, then atomically creates it without replacement.
 * @author Shinobu
 */
public final class PickSubsetMatrixGenerator{

	/** Accepts out=FILE, nodes=COUNT (at least two), and optional seed=LONG.
	 * @param args Generator options; existing outputs are never overwritten
	 * @throws IOException If opening, writing or closing the owned output fails */
	public static void main(final String[] args) throws IOException{
		String output=null;
		int nodes=0;
		long seed=1;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(key.equals("out")){output=pair[1];}
			else if(key.equals("nodes")){nodes=Integer.parseInt(pair[1]);}
			else if(key.equals("seed")){seed=Long.parseLong(pair[1]);}
			else{throw new IllegalArgumentException("Unknown matrix-generator option: "+arg);}
		}
		require(nodes>=2, "Set nodes to at least two for a complete comparison matrix.");
		require(output!=null && !output.isEmpty(), "Set out to a new raw matrix file.");
		final FileFormat ff=FileFormat.testOutput(output, FileFormat.TXT, null, false, false, false, false);
		require(ff!=null && ff.raw() && !ff.stdio(), "The matrix generator requires a raw file, not compression or standard output.");
		final Path destination=Paths.get(ff.name()).toAbsolutePath().normalize();
		require(Tools.testOutputFiles(false, false, false, destination.toString()), "Output already exists or is not writable: "+destination);
		final byte[][] names=new byte[nodes][];
		for(int i=0; i<nodes; i++){
			names[i]=("node"+i).getBytes(StandardCharsets.US_ASCII);
		}
		final Random random=new Random(seed);
		final ByteBuilder row=new ByteBuilder(256);
		long rows=0, bytes=HEADER.length;
		//CREATE_NEW protects the actual destination even if it appears after preflight.
		try(final OutputStream file=Files.newOutputStream(destination, StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
				final BufferedOutputStream stream=new BufferedOutputStream(file, 1<<16)){
			stream.write(HEADER);
			for(int query=0; query<nodes; query++){
				for(int ref=0; ref<nodes; ref++){
					final boolean self=query==ref;
					final int ani=self ? 100000 : 90000+random.nextInt(10000);
					final int kid=self ? 1000 : 100+random.nextInt(900);
					row.clear();
					row.append(names[query]).append('\t').append(names[ref]).append('\t');
					appendThousandths(row, ani);
					row.append(METADATA);
					appendThousandths(row, kid);
					row.append(SUFFIX);
					stream.write(row.array, 0, row.length());
					bytes+=row.length();
					rows++;
				}
			}
		}
		assert(rows==(long)nodes*nodes) : "The complete matrix requires one row per ordered node pair, including self rows; nodes="+nodes+", rows="+rows;
		System.err.println("Generated nodes="+nodes+", rows="+rows+", bytes="+bytes+", seed="+seed);
	}

	/** Appends a nonnegative fixed-point value with exactly three decimal digits.
	 * @param row Reused output row
	 * @param value Integer number of thousandths */
	private static void appendThousandths(final ByteBuilder row, final int value){
		assert(value>=0) : "The decimal digit formatter requires nonnegative generated ANI/KID thousandths: "+value;
		final int fraction=value%1000;
		row.append(value/1000).append('.');
		row.append((byte)('0'+fraction/100));
		row.append((byte)('0'+(fraction/10)%10));
		row.append((byte)('0'+fraction%10));
	}

	/** Rejects invalid fixture configuration independently of assertion settings. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Twelve-column comparison-table header, written once. */
	private static final byte[] HEADER=("#Query\tRef\tANI\tQSize\tRefSize\tQBases\tRBases\tQTaxID\tRTaxID\tKID\tWKID\tSSU\n").getBytes(StandardCharsets.US_ASCII);
	/** Constant metadata following ANI and preceding KID. */
	private static final byte[] METADATA="\t1\t1\t1\t1\t-1\t-1\t".getBytes(StandardCharsets.US_ASCII);
	/** Unused WKID and absent SSU, followed by the row terminator. */
	private static final byte[] SUFFIX="\t0\t.\n".getBytes(StandardCharsets.US_ASCII);
}
