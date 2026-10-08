package ml;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Explicit native weight-precision conversion, preserving biases and topology.
 * @author Brian, Yoimiya
 */
public final class CellNetWeightBits {

	/** Loads a complete native network and atomically publishes a fresh output file. */
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
		CellNet.codingA48Out=true;
		CellNet.OUT_DENSE=false; CellNet.OUT_SPARSE=false;
		final ByteBuilder bytes=net.toBytes();
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

	private CellNetWeightBits(){}
}
