package dna;

import java.util.Arrays;

import fileIO.ByteFile;

/** Immutable elongation and initiation assignments for one genetic code.
 * Codon indices use BBTools' A=0,C=1,G=2,T=3 encoding, not NCBI's TCAG order.
 * This object does not change AminoAcid's process-wide translation tables.
 * @author Keqing
 */
public final class GeneticCode {

	/** Copies validated tables; no mutable array escapes this object. */
	private GeneticCode(final int id_, final String name_, final byte[] residues_, final boolean[] starts_){
		if(residues_.length!=64 || starts_.length!=64){
			throw new IllegalArgumentException("A genetic code needs exactly 64 codon assignments");
		}
		id=id_;
		name=name_;
		residues=residues_.clone();
		starts=starts_.clone();
		int startCount=0, stopCount=0;
		for(int i=0; i<64; i++){
			if(RESIDUES.indexOf(residues[i])<0){
				throw new IllegalArgumentException("Invalid amino-acid assignment for codon index "+i);
			}
			if(starts[i]){startCount++;}
			if(residues[i]=='*'){
				stopCount++;
				if(starts[i]){throw new IllegalArgumentException("A stop codon cannot also initiate: "+AminoAcid.codonToString(i));}
			}
		}
		if(startCount==0 || stopCount==0){
			throw new IllegalArgumentException("A context-independent calling code needs starts and stops: starts="
				+startCount+", stops="+stopCount);
		}
	}

	/** NCBI table ID, or 0 for a custom table (never a claimed NCBI code). */
	public int id(){return id;}
	public String name(){return name;}

	/** Supported context-independent codes; unsupported IDs never fall back. */
	public static GeneticCode forTable(final int id){
		switch(id){
			case 11: return BACTERIAL;
			case 4: return MYCOPLASMA;
			case 25: return GRACILIBACTERIA;
			default: throw new IllegalArgumentException("Unsupported translation table "+id+"; supported: 4,11,25, or a complete custom TSV");
		}
	}

	/** Elongation residue; -1 denotes an ambiguous codon and translates to X. */
	public byte aminoAcid(final int codon){
		checkCodon(codon);
		return codon<0 ? (byte)'X' : residues[codon];
	}

	public boolean isStart(final int codon){
		checkCodon(codon);
		return codon>=0 && starts[codon];
	}

	public boolean isStop(final int codon){
		checkCodon(codon);
		return codon>=0 && residues[codon]=='*';
	}

	/** Only a known, complete initiation site receives M. Edge fragments pass false. */
	public byte translateCodon(final int codon, final boolean initiation){
		if(!initiation){return aminoAcid(codon);}
		if(!isStart(codon)){
			throw new IllegalArgumentException("Cannot initiate at non-start codon "+AminoAcid.codonToString(codon)+" under "+name);
		}
		return (byte)'M';
	}

	/** Three DNA bases in biological orientation; ambiguous/non-DNA bytes yield -1. */
	public static int codon(final byte[] bases, final int offset){
		if(bases==null || offset<0 || offset>bases.length-3){
			throw new IllegalArgumentException("A codon needs three bases at offset "+offset+"; length="+(bases==null ? -1 : bases.length));
		}
		final int a=baseNumber(bases[offset]), b=baseNumber(bases[offset+1]), c=baseNumber(bases[offset+2]);
		return (a|b|c)<0 ? -1 : (a<<4)|(b<<2)|c;
	}

	/** Small configuration-only convenience; strict length prevents truncated keys. */
	public static int codon(final String bases){
		if(bases==null || bases.length()!=3){throw new IllegalArgumentException("Codon must have exactly three DNA bases: "+bases);}
		final int a=baseNumber(bases.charAt(0)), b=baseNumber(bases.charAt(1)), c=baseNumber(bases.charAt(2));
		return (a|b|c)<0 ? -1 : (a<<4)|(b<<2)|c;
	}

	private static int baseNumber(final int base){
		// NCBI's published tables use DNA. Do not let U create a duplicate alias of T in custom tables.
		return base<0 || base>=AminoAcid.baseToNumber.length || base=='U' || base=='u' ? -1 : AminoAcid.baseToNumber[base];
	}

	private static void checkCodon(final int codon){
		if(codon< -1 || codon>63){throw new IllegalArgumentException("Codon index must be -1 or 0..63, got "+codon);}
	}

	/** Complete TSV: codon, amino_acid, start. Header required, '#' comments allowed.
	 * Exactly 64 distinct DNA codons, one uppercase residue or '*', and 0/1 start.
	 * Fixed-width data rows deliberately reject extra fields, whitespace and partial tables.
	 */
	public static GeneticCode load(final String path){
		if(path==null || path.isEmpty()){throw new IllegalArgumentException("Custom genetic-code path is required");}
		final byte[] residues=new byte[64];
		final boolean[] starts=new boolean[64], seen=new boolean[64];
		final ByteFile bf=ByteFile.makeByteFile1(path, false);
		Throwable failure=null;
		boolean header=false;
		int count=0;
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0 || line[0]=='#'){continue;}
				if(!header){
					if(!Arrays.equals(line, HEADER)){throw invalid(path, lineNumber, "expected header codon<TAB>amino_acid<TAB>start");}
					header=true;
					continue;
				}
				if(line.length!=7 || line[3]!='\t' || line[5]!='\t'){
					throw invalid(path, lineNumber, "expected a three-base codon, one residue and 0/1 start separated by tabs");
				}
				final int code=codon(line, 0);
				if(code<0){throw invalid(path, lineNumber, "codon must contain only A/C/G/T");}
				if(seen[code]){throw invalid(path, lineNumber, "duplicate codon "+AminoAcid.codonToString(code));}
				if(RESIDUES.indexOf(line[4])<0){throw invalid(path, lineNumber, "residue must be a canonical uppercase amino acid or '*'");}
				if(line[6]!='0' && line[6]!='1'){throw invalid(path, lineNumber, "start must be 0 or 1");}
				if(line[4]=='*' && line[6]=='1'){throw invalid(path, lineNumber, "a stop codon cannot also initiate");}
				seen[code]=true;
				residues[code]=line[4];
				starts[code]=(line[6]=='1');
				count++;
			}
			if(!header || count!=64){throw invalid(path, lineNumber, "a complete code needs 64 unique codons; found "+count);}
		}catch(RuntimeException | Error e){failure=e; throw e;}
		finally{
			if(bf.close()){
				final IllegalStateException e=new IllegalStateException("Error closing custom genetic-code input "+path);
				if(failure==null){throw e;}
				failure.addSuppressed(e);
			}
		}
		return new GeneticCode(0, "custom:"+path, residues, starts);
	}

	private static IllegalArgumentException invalid(final String path, final long line, final String reason){
		return new IllegalArgumentException("Invalid genetic code "+path+" at line "+line+": "+reason);
	}

	/** NCBI AAs/Starts strings are TCAG-ordered; convert once to native packed indices. */
	private static GeneticCode fromNcbi(final int id, final String name, final String amino, final String start){
		if(amino.length()!=64 || start.length()!=64){throw new IllegalArgumentException("NCBI table "+id+" must contain 64 AAs and Starts entries");}
		final String order="TCAG";
		final byte[] residues=new byte[64];
		final boolean[] starts=new boolean[64];
		for(int i=0; i<64; i++){
			final int code=(baseNumber(order.charAt(i>>4))<<4)
				|(baseNumber(order.charAt((i>>2)&3))<<2)|baseNumber(order.charAt(i&3));
			final char flag=start.charAt(i);
			if(flag!='M' && flag!='-' && flag!='*'){throw new IllegalArgumentException("Invalid NCBI initiation flag at "+i);}
			if((flag=='*')!=(amino.charAt(i)=='*')){throw new IllegalArgumentException("Context-dependent stop assignment is unsupported for NCBI table "+id);}
			residues[code]=(byte)amino.charAt(i);
			starts[code]=(flag=='M');
		}
		return new GeneticCode(id, name, residues, starts);
	}

	private final int id;
	private final String name;
	private final byte[] residues;
	private final boolean[] starts;
	private static final String RESIDUES="ACDEFGHIKLMNPQRSTVWY*";
	private static final byte[] HEADER="codon\tamino_acid\tstart".getBytes(java.nio.charset.StandardCharsets.US_ASCII);

	// NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi, tables 11/4/25, read 2026-10-01.
	// AAs are elongation; Starts 'M' permits methionine initiation and '*' marks unconditional stop.
	private static final GeneticCode BACTERIAL=fromNcbi(11, "Bacterial, Archaeal and Plant Plastid",
		"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
		"---M------**--*----M------------MMMM---------------M------------");
	private static final GeneticCode MYCOPLASMA=fromNcbi(4, "Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma",
		"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
		"--MM------**-------M------------MMMM---------------M------------");
	private static final GeneticCode GRACILIBACTERIA=fromNcbi(25, "Candidate Division SR1 and Gracilibacteria",
		"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
		"---M------**-----------------------M---------------M------------");
}
