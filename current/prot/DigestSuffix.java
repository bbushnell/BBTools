package prot;

import java.io.FileInputStream;
import java.io.InputStream;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.util.HashMap;

/** Canonical 80-bit textual suffix of a SHA-256 digest. Full digest text is never constructed. */
final class DigestSuffix {

	static final int HEX_LENGTH=20;
	static final int BYTE_LENGTH=10;
	private static final char[] HEX="0123456789abcdef".toCharArray();

	private DigestSuffix(){}

	/** Hashes the file's stored bytes, including compression bytes, and returns its lowercase sha80. */
	static String file(final String path){
		try{
			final MessageDigest md=MessageDigest.getInstance("SHA-256");
			final byte[] buffer=new byte[1<<20];
			try(InputStream in=new DigestInputStream(new FileInputStream(path), md)){
				while(in.read(buffer)>=0){}
			}
			return fromDigest(md.digest());
		}catch(Exception e){throw new RuntimeException("Failed to hash "+path, e);}
	}

	/** Returns the lowercase sha80 of the supplied byte sequence. */
	static String bytes(final byte[] bytes){
		try{return fromDigest(MessageDigest.getInstance("SHA-256").digest(bytes));}
		catch(Exception e){throw new RuntimeException("Failed to hash byte content", e);}
	}

	/** Encodes the last ten digest bytes; rejects null or shorter input. */
	static String fromDigest(final byte[] digest){
		if(digest==null || digest.length<BYTE_LENGTH){throw new IllegalArgumentException("Digest is too short");}
		final char[] out=new char[HEX_LENGTH];
		final int start=digest.length-BYTE_LENGTH;
		for(int i=0; i<BYTE_LENGTH; i++){
			final int x=digest[start+i]&0xff;
			out[2*i]=HEX[x>>>4]; out[2*i+1]=HEX[x&15];
		}
		return new String(out);
	}

	/** Validates a lowercase sha80 or legacy full digest, returning only its final twenty characters. */
	static String normalizeRecorded(final String text, final String label){
		if(text==null){throw new IllegalArgumentException("Missing digest for "+label);}
		if(text.length()!=HEX_LENGTH && text.length()!=64){
			throw new IllegalArgumentException("Digest for "+label+" must contain 20 suffix hex characters or one legacy 64-character value");
		}
		for(int i=0; i<text.length(); i++){
			final char c=text.charAt(i);
			if(!((c>='0' && c<='9') || (c>='a' && c<='f'))){
				throw new IllegalArgumentException("Digest for "+label+" is not lowercase hexadecimal");
			}
		}
		return text.substring(text.length()-HEX_LENGTH);
	}

	/** Requires exactly twenty lowercase hexadecimal characters, with no legacy full-digest fallback. */
	static String requireSuffix(final String text, final String label){
		if(text==null || text.length()!=HEX_LENGTH){
			throw new IllegalArgumentException(label+" must be exactly 20 lowercase hexadecimal characters");
		}
		return normalizeRecorded(text, label);
	}

	/** Reads a schema-appropriate pin; schema 4 values must already have been normalized by the loader. */
	static String requiredProfileHeader(final HashMap<String, String> header, final String base, final String path){
		final String schema=header.get("schema_version");
		if("5".equals(schema) || "6".equals(schema)){
			if(header.containsKey(base+"_sha256")){throw new IllegalArgumentException("Modern schema forbids legacy digest header for "+base+": "+path);}
			final String value=header.get(base+"_sha80");
			if(value==null){throw new IllegalArgumentException("Modern schema missing digest header for "+base+": "+path);}
			return requireSuffix(value, base+"_sha80");
		}
		if("4".equals(schema)){
			if(header.containsKey(base+"_sha80")){throw new IllegalArgumentException("Schema 4 forbids sha80 digest header for "+base+": "+path);}
			final String value=header.get(base+"_sha256");
			if(value==null || value.length()!=HEX_LENGTH){
				throw new IllegalArgumentException("Schema 4 legacy digest was not normalized at load for "+base+": "+path);
			}
			return requireSuffix(value, base+"_sha256 legacy suffix");
		}
		throw new IllegalArgumentException("Digest profile has unsupported schema_version: "+path);
	}

	/** Accepts exactly one modern or legacy header, validates it, and returns the normalized suffix. */
	static String requiredCompatibleHeader(final HashMap<String, String> header, final String base, final String path){
		final String modern=header.get(base+"_sha80"), legacy=header.get(base+"_sha256");
		if(modern!=null && legacy!=null){throw new IllegalArgumentException("Both modern and legacy digest headers exist for "+base+": "+path);}
		if(modern!=null){return requireSuffix(modern, base+"_sha80");}
		if(legacy!=null){return normalizeRecorded(legacy, base+"_sha256");}
		throw new IllegalArgumentException("File missing required digest header for '"+base+"': "+path);
	}

	/** Decodes a validated sha80 into ten bytes in hexadecimal display order. */
	static byte[] decodeSuffix(final String text, final String label){
		final String suffix=requireSuffix(text, label);
		final byte[] out=new byte[BYTE_LENGTH];
		for(int i=0; i<out.length; i++){
			out[i]=(byte)((Character.digit(suffix.charAt(2*i), 16)<<4)|Character.digit(suffix.charAt(2*i+1), 16));
		}
		return out;
	}
}
