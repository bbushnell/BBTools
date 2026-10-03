package prot;

import java.util.ArrayList;
import java.util.List;

import fileIO.ByteFile;
import structures.ByteBuilder;

/**
 * Loads CallGenes' {@code outa=} protein FASTA into search-ready ProteinSequence objects with
 * UNIQUE per-gene IDs, mirroring scripts/04_fix_headers_v3.sh's exact convention (the training
 * pipeline's own fix for the same problem) so IDs stay consistent between training and
 * inference.
 * <p>
 * CallGenes emits every gene on one contig with an IDENTICAL pre-tab header (the contig's own
 * name+description) -- genes on that contig are distinguished only by the tab-separated
 * strand/position fields that follow. A generic "first token is the ID" FASTA loader (e.g.
 * {@link ProteinSearch#readFasta}) therefore collapses every gene on a contig onto ONE id,
 * which {@link ProteinSearcher#search} correctly rejects as a duplicate-identifier error rather
 * than silently mis-attributing hits (real failure found 2026-08-29 running this on CallGenes'
 * actual output: 4117 genes all shared the id "tid|511145|NC_000913.3").
 * <p>
 * The fix, matching 04_fix_headers_v3.sh exactly: the unique-ification key is the FIRST
 * TAB-separated field (not the first whitespace token -- this includes the description, with
 * internal spaces turned to underscores), and the id is {@code <key>_g<N>} where N is a
 * 0-based counter that resets whenever the key changes from the previous record. This assumes
 * (as the training pipeline does, verified there on real corpora) that one contig's genes
 * appear as a CONTIGUOUS run in the file -- true for CallGenes' own per-contig-sequential
 * output, single-threaded or not.
 *
 * @author Eru
 */
public class CallGenesProteinLoader {

	/**
	 * Loads one CallGenes {@code outa=} protein FASTA, assigning unique per-gene ids.
	 * @param fname Protein FASTA path (CallGenes outa= output)
	 * @return Search-ready protein sequences, one per gene, ids matching 04_fix_headers_v3.sh's
	 * convention exactly
	 */
	public static List<ProteinSequence> load(final String fname){
		final ArrayList<ProteinSequence> list=new ArrayList<ProteinSequence>();
		final ByteFile bf=ByteFile.makeByteFile(fname, false);
		String id=null;
		final IdMaker idMaker=new IdMaker();
		int skippedStops=0;
		final ByteBuilder seq=new ByteBuilder();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='>'){
				if(id!=null && ProteinSearch.addOrSkipMalformedStop(list, id, seq.toBytes())){skippedStops++;}
				id=idMaker.next(line);
				seq.clear();
			}else{
				seq.append(line);
			}
		}
		if(id!=null && ProteinSearch.addOrSkipMalformedStop(list, id, seq.toBytes())){skippedStops++;}
		bf.close();
		if(skippedStops>0){
			System.err.println("WARNING: skipped "+skippedStops+" CallGenes record(s) with a leading/internal '*' stop marker from "+fname+".");
		}
		if(list.isEmpty()){throw new RuntimeException("No sequences found in "+fname);}
		return list;
	}

	/** Stateful implementation of the 04_fix_headers_v3.sh CallGenes id contract. */
	static final class IdMaker{
		private String previousKey;
		private int counter;

		String next(final byte[] header){
			if(header==null || header.length<2 || header[0]!='>'){
				throw new IllegalArgumentException("CallGenes header must begin with '>'");
			}
			final String key=firstTabField(header);
			if(!key.equals(previousKey)){counter=0;previousKey=key;}
			assert(counter>=0) : "CallGenes gene counter overflow for contiguous header key "+key;
			return key+"_g"+(counter++);
		}
	}

	/**
	 * Loads and VALIDATES one CallGenes {@code outa=} protein FASTA via {@link #load(String)},
	 * additionally rejecting two malformed-input shapes that {@link #load(String)}'s purely
	 * positional unique-ification cannot itself detect (added 2026-09-02, Elly's assignment --
	 * a strengthening layer; {@link #load(String)}'s own behavior/signature is unchanged for
	 * other callers):
	 * <ul>
	 * <li>An empty (zero-residue) sequence body -- two consecutive header lines with no sequence
	 * between them, or a trailing header at end of file. {@link ProteinSequence} itself does not
	 * reject a 0-length {@code enc} array, so this would otherwise silently enter the search as a
	 * phantom always-no-hit query.</li>
	 * <li>A duplicate id AFTER unique-ification. This loader's own contract (see {@link
	 * #load(String)}'s javadoc) assumes one contig's genes form a CONTIGUOUS run; the per-key
	 * counter resets on ANY key change, including a return to an EARLIER key, so a
	 * non-contiguous repeat of the same contig header silently produces a real duplicate id
	 * (e.g. contig A, contig B, contig A again -&gt; both A-runs get {@code _g0}). {@link
	 * ProteinSearcher#checkDuplicateIds} would eventually catch this, but with a generic message
	 * that does not name the actual cause.</li>
	 * </ul>
	 * @param fname Protein FASTA path (CallGenes outa= output)
	 * @return Validated, search-ready protein sequences, one per gene
	 */
	public static List<ProteinSequence> loadValidated(final String fname){
		final List<ProteinSequence> list=load(fname);
		final java.util.HashSet<String> seenIds=new java.util.HashSet<String>();
		for(final ProteinSequence p : list){
			if(p.length()==0){
				throw new RuntimeException("Malformed CallGenes input: empty (zero-residue) "+
					"sequence for id '"+p.id+"' in "+fname+" -- consecutive '>' header lines with "+
					"no sequence between them, or a trailing header at end of file.");
			}
			if(!seenIds.add(p.id)){
				throw new RuntimeException("Malformed CallGenes input: duplicate id '"+p.id+
					"' in "+fname+" after unique-ification -- this loader assumes one contig's "+
					"genes form a CONTIGUOUS run (see load()'s javadoc); a repeated/"+
					"non-contiguous contig header violates that assumption and produces a real "+
					"duplicate id, not a hypothetical one.");
			}
		}
		return list;
	}

	/** The first tab-separated field of a header line (spaces converted to underscores,
	 *  leading '&gt;' stripped) -- the unique-ification key, matching 04_fix_headers_v3.sh's
	 *  {@code name=$1; gsub(/ /,"_",name); name=substr(name,2)} exactly. */
	private static String firstTabField(final byte[] line){
		int stop=1;//skip '>'
		while(stop<line.length && line[stop]!='\t'){stop++;}
		final StringBuilder sb=new StringBuilder(stop-1);
		for(int i=1; i<stop; i++){
			final char c=(char)line[i];
			sb.append(c==' ' ? '_' : c);
		}
		return sb.toString();
	}
}
