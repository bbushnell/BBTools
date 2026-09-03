package prot;

import java.util.Arrays;

import dna.AminoAcid;

/**
 * Shared amino-acid alphabet and on-the-fly translation definition.
 *
 * <p>The class deliberately stores both the residue-to-class table and the
 * representative letters.  Consumers therefore encode raw FASTA residues
 * without first rewriting the input sequence, while emitted k-mers remain
 * human-readable and can be re-encoded by another consumer.</p>
 */
public final class ReducedAlphabet {

	private static final String AA20="A/R/N/D/C/Q/E/G/H/I/L/K/M/F/P/S/T/W/Y/V";
	private static final String LEGACY="AILMV/C/DE/FWY/GP/HKR/NQST";
	private static final String C6="AST/C/DEHKNQR/FWY/GP/ILMV";
	private static final String C7="AST/C/DEKNQR/FWY/GP/H/ILMV";
	private static final String C8="AST/C/DEKNQR/FY/GP/H/ILMV/W";
	private static final String C9="AST/C/DN/EKQR/FY/GP/H/ILMV/W";
	private static final String C12="AT/C/DN/EKQR/F/GP/H/ILV/M/S/W/Y";
	private static final String C14="AT/C/D/EN/F/G/H/ILV/KQR/M/P/S/W/Y";
	private static final String MI7="ADEHIKLMNQRSTV/C/F/G/P/W/Y";

	private ReducedAlphabet(final String name_, final String groups_, final String letters){
		name=name_; groups=groups_; symbols=letters; classes=symbols.length();
		if(classes<1){throw new IllegalArgumentException("Alphabet has no symbols: "+name);}
		if(classes>26){throw new IllegalArgumentException("Alphabet has too many symbols: "+classes);}
		map=new int[128]; Arrays.fill(map, -1);
		final String[] tokens=groups.split("/", -1);
		if(tokens.length!=classes){throw new IllegalArgumentException("Alphabet/key class count mismatch: alphabet="+
				letters+" key_groups="+groups);}
		for(int group=0; group<tokens.length; group++){
			final String token=tokens[group];
			if(token.length()==0){throw new IllegalArgumentException("Empty group in key: "+groups);}
			if(token.charAt(0)!=symbols.charAt(group)){
				throw new IllegalArgumentException("Group representative "+token.charAt(0)+
					" does not match alphabet symbol "+symbols.charAt(group)+" in "+name);
			}
			for(int i=0; i<token.length(); i++){
				final char c=Character.toUpperCase(token.charAt(i));
				if(c<'A' || c>'Z'){throw new IllegalArgumentException("Invalid residue in key: "+token.charAt(i));}
				if(map[c]>=0 && map[c]!=group){throw new IllegalArgumentException("Residue appears in multiple key groups: "+c);}
				map[c]=group; map[Character.toLowerCase(c)]=group;
			}
		}
	}

	/** Builds a named or explicit amino alphabet. */
	public static ReducedAlphabet parse(final String alphabetSpec, final String keySpec){
		String raw=(alphabetSpec==null || alphabetSpec.length()==0 ? "amino" : alphabetSpec.trim());
		if(raw.equalsIgnoreCase("nt")){throw new IllegalArgumentException("ReducedAlphabet cannot encode alphabet=nt");}
		String key=(keySpec==null ? null : keySpec.trim());
		String groups;
		String letters;
		String name;
		if(key==null || key.length()==0){
			groups=namedGroups(raw);
			if(groups!=null){letters=representatives(groups); name=canonicalName(raw);}
			else{
				letters=normalizeLetters(raw);
				groups=singletonGroups(letters); name=letters;
			}
		}else{
			groups=keyGroups(key);
			if(groups==null){throw new IllegalArgumentException("Unknown alphabet key: "+key);}
			if(raw.equalsIgnoreCase("amino")){letters=representatives(groups);}
			else if(namedGroups(raw)!=null){
				letters=representatives(groups);
				if(!letters.equals(representatives(namedGroups(raw)))){
					throw new IllegalArgumentException("alphabet and key disagree: "+raw+" / "+key);
				}
			}else{letters=normalizeLetters(raw);}
			name=raw+"/"+key;
		}
		return new ReducedAlphabet(name, groups, letters);
	}

	/** Convenience for the assay's historical names. */
	public static ReducedAlphabet named(final String raw){return parse(raw, null);}

	private static String canonicalName(final String raw){
		if(raw.equalsIgnoreCase("amino8")){return "legacy";}
		return raw.toLowerCase();
	}

	private static String namedGroups(final String raw){
		final String n=raw.toLowerCase();
		if(n.equals("aa20") || n.equals("amino")){return AA20;}
		if(n.equals("legacy") || n.equals("amino8")){return LEGACY;}
		if(n.equals("c6")){return C6;}
		if(n.equals("c7")){return C7;}
		if(n.equals("c8")){return C8;}
		if(n.equals("c9")){return C9;}
		if(n.equals("c12")){return C12;}
		if(n.equals("c14")){return C14;}
		if(n.equals("mi7")){return MI7;}
		return null;
	}

	private static String keyGroups(final String raw){
		final String named=namedGroups(raw);
		return named!=null ? named : (raw.indexOf('/')>=0 ? raw.toUpperCase() : null);
	}

	private static String normalizeLetters(final String raw){
		final String letters=raw.toUpperCase();
		if(letters.length()==0){throw new IllegalArgumentException("Alphabet is empty");}
		final boolean[] seen=new boolean[26];
		for(int i=0; i<letters.length(); i++){
			final char c=letters.charAt(i);
			if(c<'A' || c>'Z'){throw new IllegalArgumentException("Invalid alphabet symbol: "+c);}
			if(seen[c-'A']){throw new IllegalArgumentException("Duplicate alphabet symbol: "+c);}
			seen[c-'A']=true;
		}
		return letters;
	}

	private static String representatives(final String groups){
		final String[] tokens=groups.split("/", -1);
		final StringBuilder sb=new StringBuilder(tokens.length);
		for(String token : tokens){if(token.length()==0){throw new IllegalArgumentException("Empty key group");} sb.append(token.charAt(0));}
		return sb.toString();
	}

	private static String singletonGroups(final String letters){
		final StringBuilder sb=new StringBuilder(letters.length()*2);
		for(int i=0; i<letters.length(); i++){if(i>0){sb.append('/');} sb.append(letters.charAt(i));}
		return sb.toString();
	}

	public int code(final byte residue){return residue>=0 && residue<map.length ? map[residue] : -1;}
	/** Maps BBTools' encoded standard amino-acid number (0..19) on the same table. */
	public int codeEncoded(final int encoded){
		return encoded>=0 && encoded<20 ? code(AminoAcid.numberToAcid[encoded]) : -1;
	}
	public char symbol(final int code){if(code<0 || code>=symbols.length()){throw new IllegalArgumentException("Invalid alphabet code: "+code);} return symbols.charAt(code);}
	public String symbols(){return symbols;}
	/** Slash-separated groups, retained for consumers that need the canonical definition. */
	public String groups(){return groups;}
	public String name(){return name;}
	public int classes(){return classes;}
	public int bits(){return bitsFor(classes);}
	public boolean contains(final byte residue){return code(residue)>=0;}

	public static int bitsFor(final int states){
		if(states<1){throw new IllegalArgumentException("Alphabet must have at least one symbol");}
		int bits=0; long capacity=1;
		while(capacity<states){capacity<<=1; bits++;}
		return Math.max(1, bits);
	}

	private final String name;
	private final String groups;
	private final String symbols;
	private final int classes;
	private final int[] map;
}
