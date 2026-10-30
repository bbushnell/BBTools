package test.simd;

import java.util.Random;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import structures.ByteBuilder;

/** Generates deterministic uncompressed FASTQ and its independently normalized output.
 * Immutable banks contain varied headers, lengths, bases and qualities; bank IDs
 * repeat deliberately. Expected output applies IUPAC masking, U-to-T, uppercase
 * and zero quality for N. It does not call Read validation or SIMD transforms.
 * @author Shinobu */
public final class NormalizationFastqData{

	/** Generates a fixture and its expected all-flags output.
	 * @param args out=, expected=, reads=, and optional mixed=t|f */
	public static void main(final String[] args){
		String out=null, expected=null;
		long reads=0;
		boolean mixed=true;
		for(final String arg : args){
			final int split=arg.indexOf('=');
			if(split<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=arg.substring(0, split), value=arg.substring(split+1);
			if(key.equalsIgnoreCase("out")){out=value;}
			else if(key.equalsIgnoreCase("expected")){expected=value;}
			else if(key.equalsIgnoreCase("reads")){reads=Parse.parseKMG(value);}
			else if(key.equalsIgnoreCase("mixed")){mixed=Parse.parseBoolean(value);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(out==null || expected==null || out.equals(expected) || reads<1){
			throw new IllegalArgumentException("Distinct out/expected paths and positive reads are required");
		}
		assert((BANKS&(BANKS-1))==0) : "Bank selection uses a power-of-two mask";
		final byte[][] inputs=new byte[BANKS][], outputs=new byte[BANKS][];
		final Random random=new Random(41);
		for(int bank=0; bank<BANKS; bank++){
			final int length=150+(bank&15);
			final byte[] bases=new byte[length], qualities=new byte[length];
			final byte[] normalized=new byte[length], normalizedQualities=new byte[length];
			for(int i=0; i<length; i++){
				final int choice=random.nextInt(64);
				final byte b=(mixed && choice<8 ? MIXED[random.nextInt(MIXED.length)] : DNA[choice&3]);
				final byte q=(byte)(20+random.nextInt(21));
				bases[i]=b;
				qualities[i]=q;
				final byte upper=(b>='a' && b<='z' ? (byte)(b-32) : b);
				final byte n=(upper=='U' ? (byte)'T' :
						upper=='A' || upper=='C' || upper=='G' || upper=='T' ? upper : (byte)'N');
				normalized[i]=n;
				normalizedQualities[i]=(n=='N' ? 0 : q);
			}
			inputs[bank]=record(bank, bases, qualities);
			outputs[bank]=record(bank, normalized, normalizedQualities);
		}
		final ByteStreamWriter inputWriter=new ByteStreamWriter(out, false, false, false, FileFormat.FASTQ);
		final ByteStreamWriter expectedWriter=new ByteStreamWriter(expected, false, false, false, FileFormat.FASTQ);
		inputWriter.start();
		expectedWriter.start();
		long bytes=0, bases=0;
		try{
			for(long i=0; i<reads; i++){
				final int bank=(int)(i&(BANKS-1));
				inputWriter.print(inputs[bank]);
				expectedWriter.print(outputs[bank]);
				bytes+=inputs[bank].length;
				bases+=150+(bank&15);
			}
		}finally{
			final boolean inputError=inputWriter.poisonAndWait();
			final boolean expectedError=expectedWriter.poisonAndWait();
			if(inputError || expectedError){throw new IllegalStateException("FASTQ fixture write failed");}
		}
		System.out.println("NORMALIZATION_FASTQ reads="+reads+" bases="+bases+
				" bytes="+bytes+" banks="+BANKS+" mixed="+mixed+" seed=41");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Prevents construction of this fixture generator. */
	private NormalizationFastqData(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Encodes one immutable four-line FASTQ record with Phred+33 qualities.
	 * @param bank Repeated bank identifier
	 * @param bases Base sequence
	 * @param qualities Numeric quality scores
	 * @return Complete record including its final newline */
	private static byte[] record(final int bank, final byte[] bases, final byte[] qualities){
		assert(bases.length==qualities.length) : "FASTQ requires one quality per base";
		final ByteBuilder bb=new ByteBuilder(400);
		bb.append("@simd_bank_").append(bank).nl();
		bb.append(bases).nl().append('+').nl();
		for(final byte q : qualities){bb.append((byte)(q+33));}
		bb.nl();
		return bb.toBytes();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable input-bank count, kept a power of two. */
	private static final int BANKS=4096;
	/** Ordinary unambiguous DNA bases. */
	private static final byte[] DNA={'A', 'C', 'G', 'T'};
	/** Case, RNA and ambiguity symbols for a minority of mixed-record positions. */
	private static final byte[] MIXED={'a', 'c', 'g', 't', 'u', 'U', 'R', 'y', 'S', 'w',
			'K', 'm', 'B', 'd', 'H', 'v', 'N', 'n', 'X', 'x'};
}
