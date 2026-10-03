package prot;

import java.nio.charset.StandardCharsets;

/**
 * An in-memory protein sequence: an identifier plus its validated, BBTools-encoded
 * amino-acid residues. This is the unit of input to {@link ProteinSearcher}.
 *
 * <p>The identifier is the first whitespace-delimited token of a FASTA header
 * (the frozen protein-search contract's ID semantics); construction removes one
 * conventional terminal stop marker ({@code '*'}) before validating and encoding
 * the residues via {@link Blosum62#encode}. A leading or internal {@code '*'} is
 * not a terminator and is rejected loudly rather than reaching the aligner.</p>
 *
 * @author Eru
 */
public final class ProteinSequence {

	/** Sequence identifier (first whitespace-delimited header token). */
	public final String id;
	/** Encoded residues (values 0-19 or Blosum62.X_CODE). */
	public final byte[] enc;

	/**
	 * Builds a protein sequence from raw ASCII residues. One trailing {@code '*'}
	 * is a FASTA stop marker and is stripped; every other non-BLOSUM residue,
	 * including a leading or internal {@code '*'}, is rejected loudly.
	 * @param id Sequence identifier (must be non-empty, tab-free).
	 * @param raw Raw ASCII amino-acid bytes.
	 */
	public ProteinSequence(final String id, final byte[] raw){
		if(id==null || id.length()==0){
			throw new RuntimeException("Empty sequence identifier.");
		}
		if(id.indexOf('\t')>=0){
			throw new RuntimeException("Sequence identifier contains a tab: '"+id+"'.");
		}
		if(raw==null){
			throw new RuntimeException("Null residues for sequence '"+id+"'.");
		}
		final int length=(raw.length>0 && raw[raw.length-1]=='*' ? raw.length-1 : raw.length);
		if(length==0){
			throw new RuntimeException("Empty protein after stripping terminal stop marker for sequence '"+id+"'.");
		}
		this.id=id;
		if(length==raw.length){
			this.enc=Blosum62.encode(raw, id);
		}else{
			final byte[] trimmed=new byte[length];
			System.arraycopy(raw, 0, trimmed, 0, length);
			this.enc=Blosum62.encode(trimmed, id);
		}
	}

	/**
	 * Builds a protein sequence from a String of residues (convenience).
	 * @param id Sequence identifier.
	 * @param residues Amino-acid residues as a String.
	 */
	public ProteinSequence(final String id, final String residues){
		this(id, residues.getBytes(StandardCharsets.US_ASCII));
	}

	/** Length in residues. @return Residue count. */
	public final int length(){return enc.length;}

	/**
	 * True when raw FASTA residues contain a stop marker that is not the one
	 * conventional terminal marker accepted by this class. Loaders may use this
	 * narrow predicate to skip a malformed CallGenes record while preserving the
	 * strict constructor behavior for every other invalid symbol.
	 */
	static boolean hasNonterminalStopMarker(final byte[] raw){
		if(raw==null){return false;}
		final int end=(raw.length>0 && raw[raw.length-1]=='*' ? raw.length-1 : raw.length);
		for(int i=0; i<end; i++){
			if(raw[i]=='*'){return true;}
		}
		return false;
	}
}
