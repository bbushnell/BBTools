package ml;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/**
 * Converts a native CellNet image to an explicit reduced-precision or full-float
 * weight representation while preserving biases and topology. The CLI accepts
 * exact in=, out= and bits= keys; bits=32 preserves decoded weights but cannot
 * recover precision already lost in a reduced-precision input.
 * @author Brian, Yoimiya
 */
public final class CellNetWeightBits{

	/**
	 * Loads a complete native network, requires a fresh destination path, applies
	 * the requested 18/24/32-bit weight policy, enables A48 coding and clears the
	 * forced dense/sparse output overrides, retaining the parsed representation.
	 * Uses concise rows and explicit decimal input indices, restoring all six writer
	 * settings even if serialization fails. Calls must be serialized with other users
	 * of CellNet's global writer settings. Writes a sibling .partial file, checks the writer's
	 * completion result, then creates a hard link at the destination and removes the temporary name.
	 */
	public static void main(String[] args) throws IOException{
		String in=null, out=null;
		int bits=-1;
		for(String arg:args){
			if(arg.startsWith("in=") && in==null){in=arg.substring(3);}
			else if(arg.startsWith("out=") && out==null){out=arg.substring(4);}
			else if(arg.startsWith("bits=") && bits<0){bits=Integer.parseInt(arg.substring(5));}
			else{throw new IllegalArgumentException("Unknown or duplicate argument: "+arg);}
		}
		if(in==null || out==null || in.isEmpty() || out.isEmpty() || (bits!=18 && bits!=24 && bits!=32)){
			throw new IllegalArgumentException("Require in=network out=fresh_network bits=18|24|32");
		}
		final Path destination=Paths.get(out).toAbsolutePath().normalize();
		if(Files.exists(destination)){throw new IllegalArgumentException("Output already exists: "+destination);}
		final CellNet net=CellNetParser.load(in, false);
		if(bits<32 && net.tags.containsKey("affine_prefix_layers")){
			throw new IllegalArgumentException("Normalize the declared affine prefix to input headers before truncating trained edges; affine normalization must remain float32");
		}
		net.setWeightBits(bits);
		final ByteBuilder bytes=encode(net);
		Files.createDirectories(destination.getParent());
		final Path temporary=Files.createTempFile(destination.getParent(), ".weight-bits-", ".partial");
		try{
			final ByteStreamWriter writer=new ByteStreamWriter(temporary.toString(), true, false, false);
			writer.start();
			try{writer.print(bytes);}finally{if(writer.poisonAndWait()){throw new IOException("Weight-precision output failed");}}
			Files.createLink(destination, temporary);
		}finally{Files.deleteIfExists(temporary);}
		System.out.println("WEIGHT_BITS_PASS bits="+bits+" edges="+net.countEdges()+" bytes="+bytes.length());
	}

	/** Encodes a stable native layout without leaking temporary global writer settings. */
	private static ByteBuilder encode(final CellNet net){
		final boolean a48=CellNet.codingA48Out, concise=CellNet.CONCISE;
		final boolean dense=CellNet.OUT_DENSE, sparse=CellNet.OUT_SPARSE;
		final boolean hex=CellNet.OUT_HEX, delta=CellNet.OUT_DELTA_A48;
		try{
			CellNet.codingA48Out=true; CellNet.CONCISE=true;
			CellNet.OUT_DENSE=false; CellNet.OUT_SPARSE=false;
			CellNet.OUT_HEX=false; CellNet.OUT_DELTA_A48=false;
			return net.toBytes();
		}finally{
			CellNet.codingA48Out=a48; CellNet.CONCISE=concise;
			CellNet.OUT_DENSE=dense; CellNet.OUT_SPARSE=sparse;
			CellNet.OUT_HEX=hex; CellNet.OUT_DELTA_A48=delta;
		}
	}

	private CellNetWeightBits(){}
}
