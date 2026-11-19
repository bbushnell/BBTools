package prot.cluster;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import prot.Blosum62;
import prot.ProteinSequence;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Shared immutable protein storage, with original FASTA text retained only where encoding loses it.
 * Input ordinals include short records. IDs use the first header token; duplicate IDs
 * remain distinct records because assignment keys are input ordinals.
 * @author Keqing
 */
public final class ClusteringInput {

	/** A literal existing filename wins over comma syntax; lists retain their declared order. */
	public static String[] inputFiles(String specification){
		if(specification==null || specification.isEmpty()){
			throw new IllegalArgumentException("Protein input requires a filename or ordered comma-separated list");
		}
		final String[] inputs=Files.isRegularFile(Paths.get(specification)) ?
			new String[]{specification} : specification.split(",", -1);
		for(String input:inputs){
			if(input.isEmpty()){throw new IllegalArgumentException("An ordered protein input list cannot contain an empty filename");}
		}
		assert(inputs.length>0) : "A nonempty input specification must identify at least one source";
		return inputs;
	}

	/** Reads through the native FASTA streamer. The CLI must disable input normalization first. */
	public static ClusteringInput load(String input, long maximumReads, boolean preserveFasta){
		return load(new String[]{input}, maximumReads, preserveFasta, 0);
	}

	/** Loads files in declared order with global ordinals, cropping before encoding or sorting.
	 * A positive tailCrop requires that many biological residues in every record;
	 * an optional terminal '*' is preserved but is not counted as a residue.
	 * Cropping validates the retained tail, not the discarded prefix. No hash pass
	 * or intermediate FASTA is needed. The read limit, if supplied, is global.
	 */
	public static ClusteringInput load(String[] inputs, long maximumReads, boolean preserveFasta, int tailCrop){
		if(inputs==null || inputs.length==0 || tailCrop<0){
			throw new IllegalArgumentException("Protein input needs an ordered file list and nonnegative tailcrop");
		}
		final ArrayList<ProteinSequence> sequences=new ArrayList<ProteinSequence>();
		final ArrayList<String> headers=preserveFasta ? new ArrayList<String>() : null;
		final ArrayList<byte[]> overrides=preserveFasta ? new ArrayList<byte[]>() : null;
		for(String input:inputs){
			if(input==null || input.isEmpty()){throw new IllegalArgumentException("Empty protein input filename");}
			if(maximumReads>=0 && sequences.size()>=maximumReads){break;}
			final long remaining=maximumReads<0 ? -1 : maximumReads-sequences.size();
			final FileFormat format=FileFormat.testInput(input, FileFormat.FASTA, null, true, true);
			if(!format.fasta()){throw new IllegalArgumentException("Protein clustering requires FASTA input: "+input);}
			final Streamer stream=StreamerFactory.makeStreamer(format, null, true, remaining, false, true, 1);
			stream.start();
			try{
				for(ListNum<Read> list=stream.nextList(); list!=null; list=stream.nextList()){
					for(Read read : list.list){
						if(read==null || read.mate!=null || read.id==null || read.bases==null){
							throw new IllegalArgumentException("Expected an unpaired protein FASTA record at ordinal "+sequences.size());
						}
						final String id=firstToken(read.id);
						final byte[] raw=tail(read.bases, tailCrop, sequences.size());
						final ProteinSequence sequence=new ProteinSequence(id, raw);
						sequences.add(sequence);
						if(preserveFasta){
							headers.add(read.id.equals(id) ? null : read.id);
							overrides.add(matchesDecoded(sequence.enc, raw) ? null : raw);
						}
					}
				}
			}finally{
				if(ReadWrite.closeStream(stream)){throw new IllegalStateException("Protein FASTA input failed: "+input);}
			}
		}
		assert(!preserveFasta || headers.size()==sequences.size() && overrides.size()==sequences.size()) :
			"Multifile loading must preserve one header and residue override slot per global ordinal";
		return new ClusteringInput(sequences.toArray(new ProteinSequence[sequences.size()]),
			preserveFasta ? headers.toArray(new String[headers.size()]) : null,
			preserveFasta ? overrides.toArray(new byte[overrides.size()][]) : null);
	}

	/** Retains the requested suffix without converting case or counting a FASTA stop marker. */
	private static byte[] tail(byte[] raw, int length, int ordinal){
		if(length==0){return raw;}
		final int residues=raw.length-(raw.length>0 && raw[raw.length-1]=='*' ? 1 : 0);
		if(residues<length){
			throw new IllegalArgumentException("Protein at input ordinal "+ordinal+" has "+residues+
				" residues, below tailcrop="+length+"; cropping must not filter the population");
		}
		final int start=residues-length;
		assert(start>=0) : "A tail crop begins within the already checked biological residue interval";
		return start==0 ? raw : Arrays.copyOfRange(raw, start, raw.length);
	}

	private ClusteringInput(ProteinSequence[] sequences, String[] headers, byte[][] overrides){
		assert(headers==null || headers.length==sequences.length && overrides.length==sequences.length) :
			"Original FASTA metadata must retain the same ordinal as the encoded sequence store";
		this.sequences=sequences; this.headers=headers; this.overrides=overrides;
	}

	/** Appends one retained record, preserving its full original header and retained residue bytes. */
	public void appendFasta(int input, ByteBuilder out){
		if(headers==null){throw new IllegalStateException("Original FASTA preservation was not requested during loading");}
		out.append('>').append(headers[input]==null ? sequences[input].id : headers[input]).nl();
		final byte[] original=overrides[input], encoded=sequences[input].enc;
		final int length=original==null ? encoded.length : original.length;
		for(int i=0; i<length; i++){
			out.append(original==null ? decode(encoded[i]) : original[i]);
			if((i+1)%80==0){out.nl();}
		}
		if(length%80!=0){out.nl();}
	}

	public static byte[] decode(byte[] encoded){
		final byte[] raw=new byte[encoded.length];
		for(int i=0; i<encoded.length; i++){raw[i]=decode(encoded[i]);}
		return raw;
	}

	private static byte decode(byte encoded){
		assert(encoded>=0 && encoded<20 || encoded==Blosum62.X_CODE) : "Decoded output must use the validated protein alphabet";
		return encoded==Blosum62.X_CODE ? (byte)'X' : AminoAcid.numberToAcid[encoded];
	}

	private static boolean matchesDecoded(byte[] encoded, byte[] raw){
		if(encoded.length!=raw.length){return false;}
		for(int i=0; i<encoded.length; i++){if(decode(encoded[i])!=raw[i]){return false;}}
		return true;
	}

	private static String firstToken(String header){
		int end=0;
		while(end<header.length() && header.charAt(end)>32){end++;}
		if(end==0){throw new IllegalArgumentException("FASTA header must begin with a nonempty identifier");}
		return end==header.length() ? header : header.substring(0, end);
	}

	public final ProteinSequence[] sequences;
	private final String[] headers;
	private final byte[][] overrides;
}
