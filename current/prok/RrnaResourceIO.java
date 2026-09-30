package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.function.Consumer;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/** Shared resource helpers for rRNA endpoint inference and resource producers.
 * Keeps production loading independent of corpus, training and diagnostic tools.
 * Parsing, error behavior and first-maximum tie handling match those tools.
 * @author Raiden
 */
final class RrnaResourceIO {

	static ByteFile open(String path){return ByteFile.makeByteFile(FileFormat.testInput(path,FileFormat.TEXT,null,true,true));}
	static ByteStreamWriter writer(String path){
		final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));w.start();return w;
	}
	static void header(ByteFile in,String expected,String path){
		final byte[] h=in.nextLine();require(h!=null && new String(h,StandardCharsets.UTF_8).equals(expected),"Unexpected header: "+path);
	}
	static void fresh(String prefix,String... suffixes){for(String s:suffixes){require(!new File(prefix+s).exists(),"Output exists: "+prefix+s);}}
	static void finish(ByteStreamWriter... writers){
		boolean error=false;for(ByteStreamWriter w:writers){error|=w.poisonAndWait();}require(!error,"Parity output failed");
	}
	static void stream(String path,Consumer<Read> action){
		final Streamer st=StreamerFactory.makeStreamer(FileFormat.testInput(path,FileFormat.FASTA,null,true,true),null,true,-1,false,true,1);
		boolean done=false;st.start();
		try{for(ListNum<Read> batch;(batch=st.nextList())!=null;){for(Read r:batch.list){action.accept(r);}st.returnList(batch);}require(!st.errorState(),"FASTA read failed: "+path);done=true;}
		finally{if(!done){st.close();}}
	}
	static String first(String s){
		require(s!=null && !s.isEmpty(),"Nonempty name required");int n=0;while(n<s.length() && !Character.isWhitespace(s.charAt(n))){n++;}return s.substring(0,n);
	}
	static String field(String header,String key){
		assert(header!=null && key!=null):"Scalar metadata requires an intact FASTA header";
		final String prefix=key+"=";int found=-1,end=-1;
		for(int pos=header.indexOf(prefix);pos>=0;pos=header.indexOf(prefix,pos+prefix.length())){
			if(pos>0 && header.charAt(pos-1)!=';' && header.charAt(pos-1)!=' '){continue;}
			require(found<0,"Duplicate metadata: "+key);found=pos+prefix.length();end=found;
			while(end<header.length() && header.charAt(end)!=';' && !Character.isWhitespace(header.charAt(end))){end++;}
		}
		require(found>=0 && end>found,"Missing metadata: "+key);return header.substring(found,end);
	}
	static String optionalField(String header,String key){
		assert(header!=null && key!=null):"Optional legacy fields still require intact original provenance";
		final String prefix=key+"=";
		for(int pos=header.indexOf(prefix);pos>=0;pos=header.indexOf(prefix,pos+prefix.length())){
			if(pos==0 || header.charAt(pos-1)==';' || header.charAt(pos-1)==' '){return field(header,key);}
		}
		return null;
	}
	static int argmax(float[] output){
		require(output!=null && output.length==25,"Argmax requires all25 endpoint outputs");int best=0;
		for(int i=0;i<output.length;i++){require(Float.isFinite(output[i]),"Nonfinite network outputs cannot be graded as a prediction");if(output[i]>output[best]){best=i;}}return best;
	}
	static int base(byte b){switch(b){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':case 'U':case 'u':return 3;default:return -1;}}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private RrnaResourceIO(){}
}
