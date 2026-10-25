package test.synth;

import java.io.File;
import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Deterministic input generator and output summarizer for synth performance checks.
 * @author Jean
 */
public final class SynthPerformanceFixture{

	private static final int RMG_CONTIGS_DEFAULT=12000;
	private static final int RMG_CONTIG_LEN_DEFAULT=2000;
	private static final int FUSE_TAXA_DEFAULT=15000;
	private static final int FUSE_CONTIGS_PER_TAXON=20;
	private static final int FUSE_CONTIG_LEN=200;
	private static final int LINE=80;
	private static final byte LF=(byte)'\n';

	/** Writes benchmark FASTA inputs or summarizes benchmark FASTQ/FASTA outputs. */
	public static void main(String[] args) throws Exception{
		final Config c=parse(args);
		if(c.mode.equals("write")){writeInputs(c);}
		else if(c.mode.equals("fastq")){System.out.println(summarizeFastq(c));}
		else if(c.mode.equals("fasta")){System.out.println(summarizeFasta(c));}
		else{throw new IllegalArgumentException("Unknown mode "+c.mode);}
	}

	/** Parses BBTools-style key=value arguments. */
	private static Config parse(final String[] args){
		final Config c=new Config();
		for(final String arg : args){
			final String[] split=arg.split("=", 2);
			check(split.length==2, "Expected key=value argument: "+arg);
			final String key=split[0].toLowerCase(Locale.ROOT);
			final String value=split[1];
			if(key.equals("mode")){c.mode=value.toLowerCase(Locale.ROOT);}
			else if(key.equals("path")){c.path=value;}
			else if(key.equals("rmgcontigs")){c.rmgContigs=Integer.parseInt(value);}
			else if(key.equals("rmglen")){c.rmgLen=Integer.parseInt(value);}
			else if(key.equals("fusetaxa")){c.fuseTaxa=Integer.parseInt(value);}
			else if(key.equals("expectedseqs")){c.expectedSeqs=Long.parseLong(value);}
			else if(key.equals("expectedbases")){c.expectedBases=Long.parseLong(value);}
			else{throw new IllegalArgumentException("Unknown parameter "+arg);}
		}
		check(c.mode!=null, "Missing mode=write|fastq|fasta");
		check(c.path!=null, "Missing path=");
		check(c.rmgContigs>0 && c.rmgLen>0 && c.fuseTaxa>0, "Counts and lengths must be positive");
		return c;
	}

	/** Writes the reusable inputs under path/, replacing old /tmp copies directly. */
	private static void writeInputs(final Config c) throws IOException{
		final File dir=new File(c.path);
		check(dir.isDirectory() || dir.mkdirs(), "Could not create "+c.path);
		final String root=c.path.endsWith("/") ? c.path : c.path+"/";
		final Counts rmg=writeRandomReadsReference(root+"rmg-ref.fa", c.rmgContigs, c.rmgLen);
		final Counts fuse=writeFuseReference(root+"fuse-sorted.fa", c.fuseTaxa);
		System.out.println("file\tsequences\tbases");
		System.out.println("rmg-ref.fa\t"+rmg.sequences+"\t"+rmg.bases);
		System.out.println("fuse-sorted.fa\t"+fuse.sequences+"\t"+fuse.bases);
	}

	/** Writes many equal-length tax-tagged contigs for the RandomReadsMG hot path. */
	private static Counts writeRandomReadsReference(final String out,
			final int contigs, final int contigLen) throws IOException{
		long bases=0;
		final ByteStreamWriter bsw=new ByteStreamWriter(out, true, false, false);
		final ByteBuilder bb=new ByteBuilder(8192);
		bsw.start();
		try{
			for(int c=0; c<contigs; c++){
				bb.clear().append(">tid|").append(100000+c).append("|rmg_").append(c).nl();
				bases+=appendBases(bb, c, contigLen);
				bsw.print(bb);
			}
		}finally{
			check(!bsw.poisonAndWait(), "Error writing "+out);
		}
		return new Counts(contigs, bases);
	}

	/** Writes adjacent same-taxon groups for FuseByTaxa's parse/append/emit loop. */
	private static Counts writeFuseReference(final String out, final int taxa) throws IOException{
		long bases=0, sequences=0;
		final ByteStreamWriter bsw=new ByteStreamWriter(out, true, false, false);
		final ByteBuilder bb=new ByteBuilder(8192);
		bsw.start();
		try{
			for(int tax=0; tax<taxa; tax++){
				for(int c=0; c<FUSE_CONTIGS_PER_TAXON; c++){
					bb.clear().append(">tid|").append(200000+tax).append("|fuse_")
							.append(tax).under().append(c).nl();
					bases+=appendBases(bb, tax*FUSE_CONTIGS_PER_TAXON+c, FUSE_CONTIG_LEN);
					bsw.print(bb);
					sequences++;
				}
			}
		}finally{
			check(!bsw.poisonAndWait(), "Error writing "+out);
		}
		return new Counts(sequences, bases);
	}

	/** Appends a deterministic, low-compression wrapped FASTA sequence. */
	private static int appendBases(final ByteBuilder bb, final int id, final int length){
		for(int i=0; i<length; i++){
			bb.append(baseAt(id, i));
			if((i+1)%LINE==0){bb.nl();}
		}
		if(length%LINE!=0){bb.nl();}
		return length;
	}

	/** Returns a deterministic A/C/G/T base without storing a large Random stream. */
	private static char baseAt(final int id, final int pos){
		long x=0x9E3779B97F4A7C15L*(id+1L)+0xBF58476D1CE4E5B9L*(pos+1L);
		x^=(x>>>30); x*=0xBF58476D1CE4E5B9L;
		x^=(x>>>27); x*=0x94D049BB133111EBL;
		x^=(x>>>31);
		return "ACGT".charAt((int)(x&3));
	}

	/** Returns reads, bases and digest for a FASTQ output without retaining it. */
	private static String summarizeFastq(final Config c) throws NoSuchAlgorithmException{
		final DigestInput in=new DigestInput(c.path);
		long reads=0, bases=0;
		try{
			for(byte[] header; (header=in.nextLine())!=null; reads++){
				final byte[] seq=in.nextLine(), plus=in.nextLine(), qual=in.nextLine();
				check(seq!=null && plus!=null && qual!=null, "Truncated FASTQ record in "+c.path);
				check(header.length>0 && header[0]=='@', "FASTQ header missing @ in "+c.path);
				check(plus.length>0 && plus[0]=='+', "FASTQ plus line missing + in "+c.path);
				check(seq.length==qual.length, "FASTQ sequence/quality length mismatch in "+c.path);
				bases+=seq.length;
			}
		}finally{in.close();}
		check(c.expectedSeqs<0 || reads==c.expectedSeqs,
				"FASTQ read count "+reads+" != expected "+c.expectedSeqs);
		check(c.expectedBases<0 || bases==c.expectedBases,
				"FASTQ base count "+bases+" != expected "+c.expectedBases);
		return "FASTQ\t"+c.path+"\treads\t"+reads+"\tbases\t"+bases+"\tsha80\t"+in.sha80();
	}

	/** Returns sequences, bases and digest for a FASTA output without retaining it. */
	private static String summarizeFasta(final Config c) throws NoSuchAlgorithmException{
		final DigestInput in=new DigestInput(c.path);
		long sequences=0, bases=0;
		try{
			for(byte[] line; (line=in.nextLine())!=null; ){
				if(line.length>0 && line[0]=='>'){sequences++;}
				else{bases+=line.length;}
			}
		}finally{in.close();}
		check(c.expectedSeqs<0 || sequences==c.expectedSeqs,
				"FASTA sequence count "+sequences+" != expected "+c.expectedSeqs);
		check(c.expectedBases<0 || bases==c.expectedBases,
				"FASTA base count "+bases+" != expected "+c.expectedBases);
		return "FASTA\t"+c.path+"\tsequences\t"+sequences+"\tbases\t"+bases+"\tsha80\t"+in.sha80();
	}

	/** Test expectations stay fatal with assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final class Counts{
		Counts(final long sequences_, final long bases_){
			sequences=sequences_;
			bases=bases_;
		}
		final long sequences, bases;
	}

	private static final class Config{
		String mode=null, path=null;
		int rmgContigs=RMG_CONTIGS_DEFAULT, rmgLen=RMG_CONTIG_LEN_DEFAULT;
		int fuseTaxa=FUSE_TAXA_DEFAULT;
		long expectedSeqs=-1, expectedBases=-1;
	}

	/** ByteFile wrapper that updates a SHA-256 digest with LF-normalized line bytes. */
	private static final class DigestInput{
		DigestInput(final String path_) throws NoSuchAlgorithmException{
			path=path_;
			in=ByteFile.makeByteFile(path, false);
			digest=MessageDigest.getInstance("SHA-256");
		}

		byte[] nextLine(){
			final byte[] line=in.nextLine();
			if(line!=null){
				digest.update(line);
				digest.update(LF);
			}
			return line;
		}

		void close(){check(!in.close(), "Error reading "+path);}

		String sha80(){
			final byte[] bytes=digest.digest();
			final ByteBuilder bb=new ByteBuilder(20);
			for(int i=bytes.length-10; i<bytes.length; i++){
				final int b=bytes[i]&0xFF;
				if(b<16){bb.append('0');}
				bb.append(Integer.toHexString(b));
			}
			return bb.toString();
		}

		final String path;
		final ByteFile in;
		final MessageDigest digest;
	}
}
