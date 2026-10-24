package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Observes downloaded prokaryotic assemblies and validates every source GFF interval.
 * The orchestration producer must verify sources= hashes before and after execution.
 * This adapter binds those hashes; it neither computes biological truth nor repairs inputs.
 * @author Brian Bushnell, Raiden */
public final class NcrnaFreshPanelValidation {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaFreshPanelValidation.class, false);
		try{
			String selection=null, downloads=null, sources=null, out=null;int expected=200;
			for(String arg:pp.args){
				final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
				final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
				if(key.equals("selection")){selection=value;}else if(key.equals("downloads")){downloads=value;}
				else if(key.equals("sources")){sources=value;}else if(key.equals("outprefix")){out=value;}
				else if(key.equals("expected")){expected=Integer.parseInt(value);}else{throw new IllegalArgumentException("Unknown flag: "+key);}
			}
			require(selection!=null && downloads!=null && sources!=null && out!=null && expected>0,
				"Bind frozen selection, downloaded assembly pairs, verified source hashes and fresh output");
			Shared.setThreads(1);ReadWrite.setZipThreads(1);ReadWrite.ZIPLEVEL=1;
			final NcrnaFreshPanelValidation x=new NcrnaFreshPanelValidation();
			x.hashes=EukPanelGffValidation.hashes(sources);
			x.bind(selection, null);x.bind(downloads, null);
			x.loadSelection(selection, expected);x.loadDownloads(downloads);x.write(out);
			pp.outstream.println("FRESH_PANEL_VALIDATION_COMPLETE assemblies="+x.selected.size()+" no_truth_admission=true no_coordinate_repair=true");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private void bind(String path, String expected){
		assert(hashes!=null):"All observations must use the producer-verified source manifest";
		final String actual=hashes.get(path);
		require(path.startsWith("/") && actual!=null && (expected==null || actual.equals(expected)),
			"Missing or conflicting verified source binding: "+path);
	}

	private void loadSelection(String path, int expected){
		assert(selected.isEmpty()):"Load the frozen panel once, without merging independent selections";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try{
			RrnaResourceIO.header(in, SELECTION_HEADER, path);
			for(byte[] line; (line=in.nextLine())!=null; ){
				p.set(line);require(p.terms()==13, "Selection schema changed");
				final Selected s=new Selected(p.parseString(0), p.parseInt(2), p.parseString(7), p.parseLong(11));
				require(s.id.matches("GCF_[0-9]+\\.[0-9]+") && ids.add(s.id) && s.taxid>0 && s.bases>0
					&& (s.domain.equals("bacteria") || s.domain.equals("archaea"))
					&& p.termEquals("KNOWN_SOURCES_EXCLUDED_NEEDS_PROVENANCE", 12),
					"Selection must preserve unique prokaryotic assembly identity, dimensions and provenance scope");
				selected.add(s);
			}
		}finally{require(!in.close(), "Selection read failed");}
		require(selected.size()==expected, "Selected assembly count differs: "+selected.size()+" expected="+expected);
	}

	private void loadDownloads(String path){
		assert(!selected.isEmpty()):"Download indices refer to the frozen selected-row order";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);int count=0;
		try{
			RrnaResourceIO.header(in, DOWNLOAD_HEADER, path);
			for(byte[] line; (line=in.nextLine())!=null; ){
				p.set(line);require(p.terms()==7, "Download schema changed");final int index=p.parseInt(0);
				require(index>0 && index<=selected.size(), "Download index outside selected panel: "+index);
				final Selected s=selected.get(index-1);
				require(s.fasta==null && p.termEquals(s.id, 1), "Duplicate index or assembly/index disagreement: "+index);
				s.fasta=p.parseString(2);s.fastaHash=p.parseString(3);s.gff=p.parseString(4);s.gffHash=p.parseString(5);s.seal=p.parseString(6);
				require(SeedOffsetTable.validSha80(s.fastaHash) && SeedOffsetTable.validSha80(s.gffHash)
					&& SeedOffsetTable.validSha80(s.seal) && paths.add(s.fasta) && paths.add(s.gff),
					"Each downloaded assembly requires distinct FASTA/GFF paths and valid provenance hashes");
				bind(s.fasta, s.fastaHash);bind(s.gff, s.gffHash);count++;
			}
		}finally{require(!in.close(), "Download manifest read failed");}
		require(count==selected.size(), "Every selected assembly must have exactly one downloaded pair");
	}

	private void write(String prefix){
		assert(!selected.isEmpty()):"Validation must cover the entire selected panel, not a caller-selected subset";
		RrnaResourceIO.fresh(prefix, ".assemblies.tsv", ".contigs.tsv", ".features.tsv.gz", ".summary.tsv");
		final EukPanelGffValidation validator=new EukPanelGffValidation();
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".assemblies.tsv");ByteStreamWriter contigs=null;
		final ByteBuilder b=new ByteBuilder();int index=0;
		try{
			contigs=RrnaResourceIO.writer(prefix+".contigs.tsv");
			out.println("index\tassembly\tdomain\tselected_taxid\tselected_bases\tobserved_bases\tcontigs\tfasta\tfasta_sha80\tgff\tgff_sha80\tdownload_seal");
			contigs.println("domain\tassembly\tseqid\tlength\tcompartment\tnative_taxid\tselected_taxid\ttaxonomy_status\tfasta\tfasta_sha80");
			for(Selected s:selected){
				index++;final ArrayList<EukPanelGffValidation.Contig> observed=new ArrayList<EukPanelGffValidation.Contig>();
				final ObjectSet<String> literals=new ObjectSet<String>(String.class), accessions=new ObjectSet<String>(String.class);
				final int[] nativeTaxid=new int[1];final long[] bases=new long[1];
				RrnaResourceIO.stream(s.fasta, r->{
					final String id=RrnaResourceIO.first(r.id), accession=EukPanelCompartments.accession(id, nativeTaxid);
					require(r.bases!=null && r.length()>0 && literals.add(id), "Empty or duplicate FASTA contig in "+s.id+": "+id);
					require(accessions.add(accession), "Ambiguous normalized FASTA accession in "+s.id+": "+accession);
					observed.add(new EukPanelGffValidation.Contig(id, r.length(), "PROKARYOTIC", nativeTaxid[0]));bases[0]+=r.length();
				});
				require(!observed.isEmpty() && bases[0]==s.bases, "Observed FASTA bases differ from frozen catalog for "+s.id+": observed="+bases[0]+" selected="+s.bases);
				final EukPanelGffValidation.Assembly a=new EukPanelGffValidation.Assembly(s.domain, s.id, s.gff, bases[0], observed.size());
				a.gffHash=s.gffHash;
				for(EukPanelGffValidation.Contig c:observed){
					require(a.literal.put(c.id, c)==null, "Observed contig list must remain unique within "+s.id);
					EukPanelGffValidation.alias(a, EukPanelCompartments.accession(c.id, nativeTaxid), c, "VALIDATED_NATIVE_TAG_ACCESSION");
					final String taxonomy=c.taxid==0?"NOT_TAGGED":c.taxid==s.taxid?"MATCH":"MISMATCH";
					contigs.print(b.clear().append(s.domain).tab().append(s.id).tab().append(c.id).tab().append(c.length).tab().append(c.compartment).tab()
						.append(c.taxid).tab().append(s.taxid).tab().append(taxonomy).tab().append(s.fasta).tab().append(s.fastaHash).nl());
				}
				validator.assemblies.add(a);
				out.print(b.clear().append(index).tab().append(s.id).tab().append(s.domain).tab().append(s.taxid).tab().append(s.bases).tab()
					.append(bases[0]).tab().append(observed.size()).tab().append(s.fasta).tab().append(s.fastaHash).tab().append(s.gff).tab().append(s.gffHash).tab().append(s.seal).nl());
			}
		}finally{if(contigs==null){RrnaResourceIO.finish(out);}else{RrnaResourceIO.finish(out, contigs);}}
		validator.write(prefix);
	}

	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Selected{
		Selected(String i, int t, String d, long b){id=i;taxid=t;domain=d;bases=b;}
		final String id, domain;final int taxid;final long bases;
		String fasta, fastaHash, gff, gffHash, seal;
	}
	private ObjectMap<String,String> hashes;
	private final ArrayList<Selected> selected=new ArrayList<Selected>();
	public static final String SELECTION_HEADER="assembly\tpaired_assembly\ttaxid\tspecies_taxid\tcatalog_species_taxid\tphylum_taxid\tphylum\tdomain\trelease_date\tftp_directory\torganism\tbases\tprovenance_status";
	public static final String DOWNLOAD_HEADER="index\tassembly\tfasta\tfasta_sha80\tgff\tgff_sha80\tdownload_seal";
}
