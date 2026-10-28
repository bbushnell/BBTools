package prok;

import java.util.ArrayList;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Extracts unambiguous CM-supported cores from long annotated sense contexts.
 * Source annotations and every rejected or unresolved row remain in the ledger.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCoreRecovery {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaCoreRecovery.class, false);
		String audit=null, contexts=null, out=null;int expected=-1, minimum=-1, clen=-1;double ga=Double.NaN;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("audit")){audit=value;}else if(key.equals("contexts")){contexts=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("expected")){expected=Parse.parseIntKMG(value);}else if(key.equals("minsource")){minimum=Parse.parseIntKMG(value);}
			else if(key.equals("clen")){clen=Parse.parseIntKMG(value);}else if(key.equals("gabits")){ga=Double.parseDouble(value);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(audit!=null && contexts!=null && out!=null && expected>0 && minimum>0 && clen>0 && Double.isFinite(ga),
			"Bind source-attributed CM ledger, original contexts, denominator, long-source minimum and model criteria");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		final NcrnaCoreRecovery tool=new NcrnaCoreRecovery();tool.load(audit, expected, minimum, clen, ga);tool.extract(contexts, out, expected);
		Shared.closeStream(pp.outstream);
	}

	private void load(final String path, final int expected, final int minimum, final int clen, final double ga){
		assert(ordered.isEmpty()):"A core recovery consumes one fixed audit population";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, AUDIT_HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==34, "Expected the complete21source+13audit fields");
				final Locus l=new Locus(p, row);require(!byId.contains(l.id), "Duplicate source-locus identity: "+l.id);
				require(l.length>=minimum && l.stop-l.start+1==l.length && l.taxid>0 && l.start>0 && l.stop<=l.contigLength,
					"Long contiguous source geometry/taxonomy must remain valid: "+l.id);
				require(p.termEquals("ANNOTATION_FULL_COPY", 14) && p.termEquals("false", 12) && p.termEquals("false", 13)
					&& p.termEquals("false", 16) && (l.strand=='+' || l.strand=='-') && l.left>=0 && l.right>=0,
					"Recovery cohort must be source-full, nonedge and sense-oriented: "+l.id);
				if(p.termEquals("CM_FULL_SPAN_SUPPORTED", 21)){
					require(p.parseInt(23)==1 && p.parseInt(24)==1 && p.parseInt(25)==clen && p.termEquals("no", 28) && p.parseDouble(29)>=ga,
						"Accepted audit tier must bind one overlapping full-model, untruncated, GA-passing hit: "+l.id);
					l.from=p.parseInt(26);l.to=p.parseInt(27);
					require(l.from>0 && l.to>=l.from && l.to<=l.length+l.left+l.right, "CM interval exceeds original context: "+l.id);
					require(p.termEquals("EXACT", 33) || p.termEquals("HIT_WITHIN_ANNOTATION", 33) || p.termEquals("ANNOTATION_WITHIN_HIT", 33),
						"Full support requires containment geometry, never a neighboring partial overlap: "+l.id);
					l.reason=l.to-l.from+1<17?"CORE_TOO_SHORT_FOR_SEED":"CANDIDATE";
				}else{l.reason=p.parseString(21);}
				byId.put(l.id, l);ordered.add(l);
			}
		}finally{require(!in.close(), "Audit ledger read failed");}
		require(ordered.size()==expected, "Every targeted long locus must remain accounted: "+ordered.size()+" expected="+expected);
	}

	private void extract(final String path, final String prefix, final int expected){
		RrnaResourceIO.fresh(prefix, ".cores.fa", ".ledger.tsv", ".summary.tsv");
		final ByteStreamWriter fasta=RrnaResourceIO.writer(prefix+".cores.fa"), ledger=RrnaResourceIO.writer(prefix+".ledger.tsv");
		final ByteBuilder b=new ByteBuilder();final int[] counts={0, 0};
		try{
			RrnaResourceIO.stream(path, r->{
				final Locus l=byId.get(RrnaResourceIO.first(r.id));require(l!=null && !l.seen, "Foreign/duplicate context: "+r.id);
				require(RrnaResourceIO.field(r.id, "source").equals(l.genome) && RrnaResourceIO.field(r.id, "contig").equals(l.seqid)
					&& RrnaResourceIO.field(r.id, "strand").equals(String.valueOf(l.strand))
					&& Integer.parseInt(RrnaResourceIO.field(r.id, "start"))==l.start && Integer.parseInt(RrnaResourceIO.field(r.id, "stop"))==l.stop
					&& Integer.parseInt(RrnaResourceIO.field(r.id, "lflank"))==l.left && Integer.parseInt(RrnaResourceIO.field(r.id, "rflank"))==l.right
					&& r.bases!=null && r.length()==l.length+l.left+l.right, "Context provenance/length changed: "+l.id);
				l.seen=true;counts[0]++;
				if(!l.reason.equals("CANDIDATE")){return;}
				for(int i=l.from-1; i<l.to; i++){
					if(AminoAcid.baseToNumber[r.bases[i]&255]<0){l.reason="AMBIGUOUS_CORE";return;}
				}
				l.coreStart=l.strand=='+'?l.start-l.left+l.from-1:l.stop+l.left-l.to+1;
				l.coreStop=l.strand=='+'?l.start-l.left+l.to-1:l.stop+l.left-l.from+1;
				require(l.coreStart>=1 && l.coreStop<=l.contigLength && l.coreStop-l.coreStart==l.to-l.from,
					"Context-contained CM coordinates must preserve length and remain inside the real source contig: "+l.id);
				l.reason="RECOVERED_CM_CORE";counts[1]++;
				b.clear().append('>').append(l.id).append(" source=").append(l.genome).append(" tid_").append(l.taxid)
					.append(" contig=").append(l.seqid).append(" strand=").append(l.strand).append(" start=").append(l.coreStart)
					.append(" stop=").append(l.coreStop).append(" historical_start=").append(l.start).append(" historical_stop=").append(l.stop)
					.append(" recovery=CM_FULL_SPAN_SUPPORTED lflank=0 rflank=0\n");
				for(int i=l.from-1; i<l.to; i++){b.append((char)Character.toUpperCase((char)r.bases[i]));}fasta.print(b.nl());
			});
			require(counts[0]==expected, "Missing original contexts cannot silently disappear from core recovery");
			ledger.println(AUDIT_HEADER+"\trecovery\tcore_start1\tcore_stop1\tcore_length");
			for(Locus l:ordered){
				require(l.seen && !l.reason.equals("CANDIDATE"), "Every source row needs a final recovery disposition: "+l.id);
				b.clear().append(l.row).tab().append(l.reason).tab();
				if(l.reason.equals("RECOVERED_CM_CORE")){b.append(l.coreStart).tab().append(l.coreStop).tab().append(l.coreStop-l.coreStart+1);}
				else{b.append("NA\tNA\tNA");}ledger.print(b.nl());
			}
		}finally{RrnaResourceIO.finish(fasta, ledger);}
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		try{summary.println("targeted\trecovered\tunresolved_or_rejected");summary.print(b.clear().append(expected).tab().append(counts[1]).tab().append(expected-counts[1]).nl());}
		finally{RrnaResourceIO.finish(summary);}
		System.err.println("NCRNA_CORE_RECOVERY_PASS targeted="+expected+" recovered="+counts[1]+"; original annotations unchanged");
	}

	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Locus {
		Locus(final LineParser1 p, final byte[] r){
			row=r;id=p.parseString(0);genome=p.parseString(1);taxid=p.parseInt(3);seqid=p.parseString(6);
			final String strandText=p.parseString(7), interval=p.parseString(8);
			require(strandText.length()==1 && interval.indexOf(',')<0, "Recovery requires one strand and one contiguous interval");strand=strandText.charAt(0);
			final int dash=interval.indexOf('-');require(dash>0, "Missing source interval delimiter");
			start=Integer.parseInt(interval.substring(0, dash));stop=Integer.parseInt(interval.substring(dash+1));
			length=p.parseInt(9);contigLength=p.parseInt(15);left=p.parseInt(18);right=p.parseInt(19);
		}
		final String id, genome, seqid;final byte[] row;final char strand;final int taxid, start, stop, length, contigLength, left, right;
		String reason;int from, to, coreStart, coreStop;boolean seen;
	}
	public static final String AUDIT_HEADER=NcrnaCmLocusAudit.HEADER+"\tcm_state\traw_hits\toverlapping_sense_hits\tcm_model_from\tcm_model_to\tcm_seq_from\tcm_seq_to\tcm_trunc\tcm_bits\toverlap_bp\tleft_excess\tright_excess\tgeometry";
	private final ObjectMap<String,Locus> byId=new ObjectMap<String,Locus>(String.class, Locus.class);
	private final ArrayList<Locus> ordered=new ArrayList<Locus>();
}
