package prot;

import java.io.IOException;
import java.io.InputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

import fileIO.ReadWrite;
import structures.ByteBuilder;

/**
 * Exact decoded-byte access for compressed MAG-QC text resources. Embedded
 * scientific pins identify the uncompressed TSV/FASTA bytes, independently of
 * BGZF block layout or compression level. The release manifest separately pins
 * stored download bytes. NN/HBM transports and input-genome hashes do not use
 * this helper. Parsers continue to read lines through ByteFile.
 * @author Collei
 */
final class MagQCTextResource {

	/** Returns the logical content pin without normalizing whitespace or line endings. */
	static String sha80(final String path){return DigestSuffix.fromDigest(digest(path));}

	/** Full digest bytes are retained only for legacy in-memory comparison contracts. */
	static byte[] digest(final String path){
		try{
			final MessageDigest md=MessageDigest.getInstance("SHA-256");
			final byte[] buffer=new byte[1<<16];
			try(InputStream in=ReadWrite.getInputStream(path, true, false)){
				for(int len=in.read(buffer); len>=0; len=in.read(buffer)){
					if(len>0){md.update(buffer, 0, len);}
				}
			}
			return md.digest();
		}catch(IOException | NoSuchAlgorithmException e){throw new IllegalArgumentException("Cannot hash decoded text resource: "+path, e);}
	}

	/**
	 * Reads exact decoded bytes for small canonical-roster checks. This preserves
	 * the existing rejection of CRLF/missing final LF, which line parsing cannot
	 * distinguish. Large tables use streaming ByteFile parsers, not this method.
	 */
	static byte[] bytes(final String path) throws IOException{
		final ByteBuilder out=new ByteBuilder(8192);
		final byte[] buffer=new byte[1<<16];
		try(InputStream in=ReadWrite.getInputStream(path, true, false)){
			for(int len=in.read(buffer); len>=0; len=in.read(buffer)){
				if(len>Integer.MAX_VALUE-8-out.length()){throw new IOException("Text resource exceeds Java byte-array capacity: "+path);}
				if(len>0){out.append(buffer, 0, len);}
			}
		}
		return out.toBytes();
	}
	private MagQCTextResource(){}
}
