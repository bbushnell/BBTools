package prok;

import java.util.ArrayList;
import java.util.Collections;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import structures.ByteBuilder;

/** Streaming multiset comparison of C-locale sorted, header-free GFF rows.
 * All columns participate in identity; changed scores/models remain visible.
 * Does not substitute call-set equality for independent biological truth.
 * @author Raiden, Ganyu
 */
public final class GffCallSetDelta {
	public static void main(String[] args){
		String before=null,after=null,out=null,mode="delta",changes=null,loci=null,ungraded=null,audit=null,panel=null;
		for(String arg:args){final int i=arg.indexOf('=');require(i>0,"Expected key=value");
			final String k=arg.substring(0,i).toLowerCase(java.util.Locale.ROOT),v=arg.substring(i+1);
			if(k.equals("before")){before=v;}else if(k.equals("after")){after=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("mode")){mode=v;}else if(k.equals("changes")){changes=v;}else if(k.equals("loci")){loci=v;}
			else if(k.equals("ungraded")){ungraded=v;}else if(k.equals("audit")){audit=v;}else if(k.equals("panel")){panel=v;}
			else{throw new IllegalArgumentException(k);}
		}
		require(out!=null,"Output prefix required");
		if(mode.equals("bind")){
			require(before==null && after==null && changes!=null,"Binding annotates saved deltas only; do not supply caller GFFs");
			bind(changes,loci,ungraded,audit,panel,out);
		}else{
			require(mode.equals("delta") && changes==null && loci==null && ungraded==null && audit==null && panel==null,
				"Metadata input requires explicit mode=bind");
			require(before!=null && after!=null,"Two sorted GFF streams required");new GffCallSetDelta().run(before,after,out);
		}
	}
	public void run(String before,String after,String prefix){
		for(String name:new String[]{"CDS","tRNA","legacy16S","legacy23S","legacy5S","legacy18S","legacy28S","euk5S"}){get(name);}
		RrnaCmWindowParity.fresh(prefix,".summary.tsv",".changes.tsv.gz");final Cursor a=new Cursor(before),z=new Cursor(after);final ByteStreamWriter changes=RrnaSourceTargetPrep.writer(prefix+".changes.tsv.gz");final ByteBuilder b=new ByteBuilder();
		try{changes.println("change\tseqid\tsource\ttype\tstart\tstop\tscore\tstrand\tphase\tattributes");
			while(a.row!=null || z.row!=null){final int compare=a.row==null?1:z.row==null?-1:compare(a.row,z.row);
				if(compare==0){final Counts c=get(family(a.row));c.before++;c.after++;c.shared++;a.advance();z.advance();}
				else if(compare<0){final Counts c=get(family(a.row));c.before++;c.removed++;changes.print(b.clear().append("REMOVED\t").append(a.row).nl());a.advance();}
				else{final Counts c=get(family(z.row));c.after++;c.added++;changes.print(b.clear().append("ADDED\t").append(z.row).nl());z.advance();}
			}
		}finally{a.close();z.close();RrnaCmWindowParity.finish(changes);}
		Collections.sort(order);final ByteStreamWriter summary=RrnaSourceTargetPrep.writer(prefix+".summary.tsv");try{summary.println("family\tbefore\tafter\tunchanged\tremoved\tadded\texactMatch");for(String name:order){final Counts c=counts.get(name);require(c.before==c.shared+c.removed && c.after==c.shared+c.added,"Multiset comparison must conserve duplicate call multiplicity");summary.print(b.clear().append(name).tab().append(c.before).tab().append(c.after).tab().append(c.shared).tab().append(c.removed).tab().append(c.added).tab().append(c.removed==0 && c.added==0).nl());}}finally{RrnaCmWindowParity.finish(summary);}
		System.err.println("GFF_CALLSET_DELTA_PASS families="+order.size());
	}
	/** Adds source ownership from native saved grading rows only. Metadata is
	 * restricted to the delta's arm and panel; ambiguous and missing owners remain
	 * visible. No genome read, feature rematching, or biological interpretation. */
	public static void bind(String changes,String loci,String ungraded,String audit,String panel,String prefix){
		final boolean annotation=loci!=null && ungraded!=null && audit==null && panel==null;
		final boolean frozen=loci==null && ungraded==null && audit!=null
			&& ("P1".equals(panel) || "P2".equals(panel) || "P3".equals(panel));
		require(annotation || frozen,"Bind either panel loci+ungraded or frozen audit+explicit panel");
		final ObjectMap<String,Owners> wanted=new ObjectMap<String,Owners>(String.class,Owners.class);
		final LineParser1 p=new LineParser1('\t');final ByteFile first=RrnaResourceIO.open(changes);long expected=0;
		try{
			RrnaResourceIO.header(first,CHANGE_HEADER,changes);
			for(byte[] row;(row=first.nextLine())!=null;){p.set(row);checkChange(p);expected++;
				final String key=changeKey(p);if(!wanted.contains(key)){wanted.put(key,new Owners());}
			}
		}finally{require(!first.close(),"Saved delta read failed");}
		if(annotation){loadOwners(loci,0,null,wanted);loadOwners(ungraded,1,null,wanted);}
		else{loadOwners(audit,2,panel,wanted);}
		RrnaResourceIO.fresh(prefix,".bindings.tsv.gz",".binding_summary.tsv");
		final ByteFile in=RrnaResourceIO.open(changes);
		final ByteStreamWriter detail=RrnaResourceIO.writer(prefix+".bindings.tsv.gz"),summary=RrnaResourceIO.writer(prefix+".binding_summary.tsv");
		final ByteBuilder b=new ByteBuilder();final long[] counts=new long[3];long rows=0;
		try{
			RrnaResourceIO.header(in,CHANGE_HEADER,changes);detail.println("binding_status\towner_count\towners\tevidence_scope\t"+CHANGE_HEADER);
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);checkChange(p);
				final Owners owners=wanted.get(changeKey(p));require(owners!=null,"Delta changed during metadata binding");
				final int size=owners.names.size(),state=size==0?2:size==1?0:1;counts[state]++;rows++;
				final ArrayList<String> names=new ArrayList<String>();for(String name:owners.names.keys()){if(name!=null){names.add(name);}}Collections.sort(names);
				b.clear().append(BINDING_STATES[state]).tab().append(size).tab();
				if(size==0){b.append("NA");}else{for(int i=0;i<size;i++){if(i>0){b.append(',');}b.append(names.get(i));}}
				detail.print(b.tab().append(annotation?"PANEL_ANNOTATION_NATIVE_ROWS":"FROZEN_AUDIT_"+panel).tab().append(row).nl());
			}
			require(rows==expected,"Binding must preserve every changed row, including duplicates");
			summary.println("binding_status\tchanged_rows");
			for(int i=0;i<counts.length;i++){summary.print(b.clear().append(BINDING_STATES[i]).tab().append(counts[i]).nl());}
		}finally{require(!in.close(),"Saved delta reread failed");RrnaResourceIO.finish(detail,summary);}
		System.err.println("GFF_DELTA_BINDING_PASS rows="+rows+" UNIQUE="+counts[0]+" AMBIGUOUS="+counts[1]+" UNRESOLVED="+counts[2]+" original_rows_preserved=true no_genome_scan=true");
	}
	private static void loadOwners(String path,int schema,String panel,ObjectMap<String,Owners> wanted){
		assert(schema>=0 && schema<=2):"Only the three explicitly supported native grading schemas supply ownership";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final String header=schema==0?PANEL_LOCI_HEADER:schema==1?PANEL_UNGRADED_HEADER:NcrnaEndpointAudit.HEADER;
		final int columns=schema==1?13:18;
		try{
			RrnaResourceIO.header(in,header,path);
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==columns,"Native ownership schema changed");
				if(!(p.termEquals("before",0) || p.termEquals("after",0))){continue;}
				if(schema==2 && !p.termEquals(panel,1)){continue;}
				final String arm=p.parseString(0),contig=p.parseString(schema==2?3:4);
				final Owners owners=wanted.get(arm+'\t'+contig);if(owners==null){continue;}
				final String assembly=p.parseString(schema==2?2:1);
				require(!assembly.isEmpty() && !assembly.equals("NA") && assembly.indexOf(',')<0,"Native owner must be an explicit unambiguous identifier token");
				owners.names.add(assembly);
			}
		}finally{require(!in.close(),"Native ownership table read failed");}
	}
	private static void checkChange(LineParser1 p){
		require(p.terms()==10 && (p.termEquals("REMOVED",0) || p.termEquals("ADDED",0)),"Expected unchanged native ten-column GFF delta schema");
		require(!p.parseString(1).isEmpty(),"Changed row needs its original seqid");
	}
	private static String changeKey(LineParser1 p){return (p.termEquals("REMOVED",0)?"before":"after")+'\t'+p.parseString(1);}
	private static final class Owners{final ObjectSet<String> names=new ObjectSet<String>(String.class);}
	public static final String CHANGE_HEADER="change\tseqid\tsource\ttype\tstart\tstop\tscore\tstrand\tphase\tattributes";
	public static final String PANEL_LOCI_HEADER="arm\tgenome\tfamily\trole\tseqid\tstart\tstop\tstrand\tphase\tsource_line\toutcome\tpartner_source_line\tpartner_start\tpartner_stop\tpartner_phase\td5\td3\tcircular_length_modular_ends_frame";
	public static final String PANEL_UNGRADED_HEADER="arm\tgenome\tscope\tsource_line\toriginal_seqid\toriginal_source\toriginal_type\toriginal_start\toriginal_stop\toriginal_score\toriginal_strand\toriginal_phase\toriginal_attributes";
	private static final String[] BINDING_STATES={"UNIQUE","AMBIGUOUS","UNRESOLVED"};
	String family(byte[] row){parser.set(row);require(parser.terms()==9,"Expected complete GFF9 rows without headers");if(parser.termEquals("CDS",2)){return "CDS";}if(parser.termEquals("tRNA",2)){return "tRNA";}
		if(parser.termEquals("internal_transcribed_spacer_region",2)){
			// RrnaItsDeriver.Derived.appendGff writes a product token, not a family tag.
			for(int i=0;i<SPACERS.length;i++){if(parser.termStartsWith(SPACER_PREFIXES[i],8)){return SPACERS[i];}}
			throw new IllegalArgumentException("Derived ITS rows require the caller's explicit ITS/ITS1/ITS2 product");
		}
		final String family=RrnaSeedBenchmarkGrade.family(row,parser,false);if(family!=null){return family;}
		for(String name:LEGACY){if(parser.termEquals("rRNA",2) && parser.termStartsWith(name+",",8)){return "legacy"+name;}}
		throw new IllegalArgumentException("Unclassified GFF family must not disappear from non-target comparison");
	}
	Counts get(String family){Counts c=counts.get(family);if(c==null){c=new Counts();counts.put(family,c);order.add(family);}return c;}
	static int compare(byte[] a,byte[] b){for(int i=0;i<Math.min(a.length,b.length);i++){final int c=(a[i]&255)-(b[i]&255);if(c!=0){return c;}}return Integer.compare(a.length,b.length);}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Cursor{Cursor(String path){in=RrnaSourceTargetPrep.open(path);advance();}void advance(){final byte[] previous=row;row=in.nextLine();require(row==null || row.length>0 && row[0]!='#',"Input must be header-free nonempty GFF rows");require(previous==null || row==null || compare(previous,row)<=0,"Input must be C-locale byte-sorted; never silently reorder a streamed multiset");}void close(){require(!in.close(),"GFF input read failed");}final ByteFile in;byte[] row;}
	public static final class Counts{public long before,after,shared,added,removed;}
	final LineParser1 parser=new LineParser1('\t');public final ObjectMap<String,Counts> counts=new ObjectMap<String,Counts>(String.class,Counts.class);final ArrayList<String> order=new ArrayList<String>();static final String[] LEGACY={"16S","23S","5S","18S","28S"};
	static final String[] SPACERS={"ITS","ITS1","ITS2"};
	static final String[] SPACER_PREFIXES={"product:ITS,","product:ITS1,","product:ITS2,"};
}
