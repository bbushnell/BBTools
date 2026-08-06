package prot;

/**
 * An in-memory protein sequence: an identifier plus its validated, BBTools-encoded
 * amino-acid residues. This is the unit of input to {@link ProteinSearcher}.
 *
 * <p>The identifier is the first whitespace-delimited token of a FASTA header
 * (the frozen protein-search contract's ID semantics); construction validates
 * and encodes the residues via {@link Blosum62#encode}, crashing loud on any
 * illegal residue.</p>
 *
 * @author Eru
 */
public final class ProteinSequence {

	/** Sequence identifier (first whitespace-delimited header token). */
	public final String id;
	/** Encoded residues (values 0-19 or Blosum62.X_CODE). */
	public final byte[] enc;

	/**
	 * Builds a protein sequence from raw ASCII residues, validating + encoding.
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
		this.id=id;
		this.enc=Blosum62.encode(raw, id);
	}

	/**
	 * Builds a protein sequence from a String of residues (convenience).
	 * @param id Sequence identifier.
	 * @param residues Amino-acid residues as a String.
	 */
	public ProteinSequence(final String id, final String residues){
		this(id, residues.getBytes());
	}

	/** Length in residues. @return Residue count. */
	public final int length(){return enc.length;}
}
