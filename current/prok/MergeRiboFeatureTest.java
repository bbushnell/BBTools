package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import fileIO.ByteStreamWriter;
import shared.Timer;
import stream.Read;
import tax.TaxTree;

/** Black-box fixtures for the complete, reconciled MergeRibo interface.
 * @author Raiden
 */
public final class MergeRiboFeatureTest {
	public static void main(String[] args){
		check(args.length==2,"Expected create|check and fixture/output directory");
		if(args[0].equals("create")){create(args[1]);}
		else{check(args[0].equals("check"),"Unknown fixture action");verify(args[1]);}
	}
	static void create(String dir){
		check(new File(dir).mkdirs(),"Feature fixture directory must be fresh");
		final byte[][] seq=MergeRiboSelectionTest.sequences();
		write(dir+"/refs.fa",new String[]{"first_wrong","second_exact"},seq[0],seq[1]);
		write(dir+"/first.fa",new String[]{"first_wrong"},seq[0]);
		write(dir+"/exact.fa",new String[]{"exact"},seq[1]);
		write(dir+"/one.fa",new String[]{"tid|101|exact"},seq[1]);
		write(dir+"/headers.fa",new String[]{"tid|101|pipe","ncbi|102|pipe","gene source=tid_103_locus",
			"gene metadata ncbi_104","plastid_noise tid_105","no_taxid","tid_0","tid_-3"},
			seq[1],seq[1],seq[1],seq[1],seq[1],seq[1],seq[1],seq[1]);
		final int[] lengths={49,50,79,80,120,200,201,300,301,6000,6001};
		final String[] lengthNames=new String[lengths.length];final byte[][] lengthSeq=new byte[lengths.length][];
		for(int i=0;i<lengths.length;i++){
			lengthNames[i]="tid|"+(1000+i)+"|length"+lengths[i];lengthSeq[i]=repeat(seq[1],lengths[i]);
		}
		write(dir+"/lengths.fa",lengthNames,lengthSeq);
		final byte[] oneN=seq[1].clone(),twoN=seq[1].clone();oneN[10]='N';twoN[10]=twoN[11]='N';
		write(dir+"/ns.fa",new String[]{"tid|41|clean","tid|42|oneN","tid|43|twoN"},seq[1],oneN,twoN);
		write(dir+"/identity.fa",new String[]{"tid|51|wrong","tid|52|exact"},seq[0],seq[1]);
		write(dir+"/primary.fa",new String[]{"tid|61|primary","tid|62|filtered_primary"},seq[1],seq[0]);
		write(dir+"/extra.fa",new String[]{"tid|62|second_primary"},seq[1]);
		write(dir+"/alt.fa",new String[]{"tid|61|ignored_alt","tid|62|fallback","tid|63|alt_only"},seq[1],seq[1],seq[1]);
		final byte[] shortSeq=Arrays.copyOf(seq[1],80);
		write(dir+"/rank.fa",new String[]{"tid|71|short","tid|71|long"},shortSeq,seq[1]);
		write(dir+"/short_first.fa",new String[]{"short","long"},shortSeq,seq[1]);
		write(dir+"/long_first.fa",new String[]{"long","short"},seq[1],shortSeq);
		finish(writer(dir+"/empty.fa"));
		for(String type:new String[]{"16S","18S","23S","5S"}){
			final Read[] refs=ProkObject.loadConsensusSequenceType(type,true,true);
			final int count=type.equals("5S")?refs.length:1;
			final String[] names=new String[count];final byte[][] bases=new byte[count][];
			for(int i=0;i<count;i++){names[i]="tid|"+(301+i)+"|"+type;bases[i]=refs[i].bases;}
			write(dir+"/default_"+type+".fa",names,bases);
		}
		final ArrayList<byte[]> its=new ArrayList<byte[]>();
		for(String type:new String[]{"ITS_fungi","ITS_plant","ITS_animal","ITS_other"}){
			final Read[] refs=ProkObject.loadConsensusSequenceType(type,false,false);if(refs!=null && refs.length>0){its.add(refs[0].bases);}
		}
		check(!its.isEmpty(),"This default-ITS fixture requires at least one shipped ITS lineage reference");
		final String[] itsNames=new String[its.size()];for(int i=0;i<itsNames.length;i++){itsNames[i]="tid|"+(401+i)+"|ITS";}
		write(dir+"/default_ITS.fa",itsNames,its.toArray(new byte[0][]));
		writeTree(dir);
		write(dir+"/taxa.fa",new String[]{"tid|111|species_a","tid|211|strain_a","tid|112|species_b"},seq[1],seq[1],seq[1]);
		System.out.println("MERGERIBO_FEATURE_FIXTURES_CREATED");
	}
	static void writeTree(String dir){
		final ByteStreamWriter names=writer(dir+"/names.dmp"),nodes=writer(dir+"/nodes.dmp");
		final int[] ids={1,2,2157,110,111,112,211},parents={1,1,1,2,110,110,111};
		final String[] labels={"root","Bacteria","Archaea","Fixture genus","Fixture genus alpha","Fixture genus beta","Fixture strain"};
		final String[] ranks={"no rank","superkingdom","superkingdom","genus","species","species","no rank"};
		try{for(int i=0;i<ids.length;i++){
			names.println(ids[i]+"\t|\t"+labels[i]+"\t|\t\t|\tscientific name\t|");
			nodes.println(ids[i]+"\t|\t"+parents[i]+"\t|\t"+ranks[i]+"\t|\t0");
		}}finally{finish(names,nodes);}
		final TaxTree tree=TaxTree.loadTaxTree(null,dir+"/names.dmp",dir+"/nodes.dmp",null,System.err,false,false);
		TaxTree.writeTaxTree(tree,dir+"/tree.taxtree.gz",false);
	}
	static void verify(String dir){
		for(String mode:new String[]{"16S","18S","ITS","LSU","5S","5.8S","23S","25S","26S","28S","58S",
			"process16S","process18S","processITS","processLSU","process5S","process5.8S"}){
			names(dir+"/ref_"+mode+".fa","tid|101|exact");
		}
		count(dir+"/first_only.fa",0);
		for(String mode:new String[]{"16S","18S","ITS","LSU","5S","5.8S"}){
			MergeRiboSelectionTest.group(load(dir+"/votes_"+mode+".fa"),11,MergeRiboSelectionTest.sequences()[0]);
		}
		names(dir+"/headers.fa","tid|101|pipe","ncbi|102|pipe","gene source=tid_103_locus","gene metadata ncbi_104","plastid_noise tid_105");
		lengths(dir+"/bounds_5S.fa",80,120,200);
		lengths(dir+"/bounds_5.8S.fa",50,79,80,120,200,201,300);
		lengths(dir+"/bounds_LSU.fa",49,50,79,80,120,200,201,300,301,6000);
		for(String mode:new String[]{"16S","18S","ITS"}){lengths(dir+"/bounds_"+mode+".fa",49,50,79,80,120,200,201,300,301);}
		lengths(dir+"/override_before.fa",79,80,120,200,201);lengths(dir+"/override_after.fa",79,80,120,200,201);
		lengths(dir+"/last_mode.fa",80,120,200);
		names(dir+"/ns.fa","tid|41|clean","tid|42|oneN");names(dir+"/identity.fa","tid|52|exact");
		names(dir+"/alt.fa","tid|61|primary","tid|62|fallback","tid|63|alt_only");
		names(dir+"/multiple_input.fa","tid|61|primary","tid|62|second_primary","tid|63|alt_only");
		names(dir+"/short_first.fa","tid|71|short");names(dir+"/long_first.fa","tid|71|long");
		for(String type:new String[]{"16S","18S","23S","5S","ITS"}){
			count(dir+"/default_"+type+".fa",load(dir+"/fixtures/default_"+type+".fa").size());
		}
		count(dir+"/species.fa",2);count(dir+"/genus.fa",1);count(dir+"/tree_true.fa",2);count(dir+"/usetree.fa",2);
		names(dir+"/reads.fa","tid|101|pipe");names(dir+"/space path/output file.fa","tid|101|exact");
		final ArrayList<Read> dada=load(dir+"/dada_implicit.fa");check(dada.size()==2,"Implicit tree loading must preserve two species");
		for(Read r:dada){check(r.id.contains("g__Fixture genus;") && r.id.contains("s__Fixture genus "),"Expected genus/species DADA2 lineage");}
		priorError(dir+"/fixtures");
		System.out.println("MERGERIBO_ALL_FEATURES_PASS refs17_aliases bounds6 overrides headers5 filters alt multiple_inputs first_reference_length defaults5 taxonomy error_propagation");
	}
	static void priorError(String fixtures){
		final MergeRibo tool=new MergeRibo(new String[]{"in="+fixtures+"/one.fa","5S=t","ref="+fixtures+"/exact.fa","tree=f","t=1"});
		tool.errorState=true;boolean failed=false;
		try{tool.process(new Timer());}catch(RuntimeException e){failed=e.getMessage().contains("terminated in an error state");}
		check(failed,"Successful input/selection phases must not erase an earlier error");
	}
	static void names(String path,String... expected){
		final ArrayList<Read> reads=load(path);final String[] actual=new String[reads.size()];
		for(int i=0;i<actual.length;i++){actual[i]=reads.get(i).id;}
		Arrays.sort(actual);Arrays.sort(expected);check(Arrays.equals(actual,expected),"Unexpected headers in "+path+": "+Arrays.toString(actual));
	}
	static void lengths(String path,int... expected){
		final ArrayList<Read> reads=load(path);final int[] actual=new int[reads.size()];
		for(int i=0;i<actual.length;i++){actual[i]=reads.get(i).length();}
		Arrays.sort(actual);Arrays.sort(expected);check(Arrays.equals(actual,expected),"Unexpected length-filter output in "+path+": "+Arrays.toString(actual));
	}
	static void count(String path,int expected){check(load(path).size()==expected,"Wrong record count in "+path+", expected "+expected);}
	static byte[] repeat(byte[] source,int length){
		check(source.length>0 && length>0,"Fixture lengths must be positive");final byte[] out=new byte[length];
		for(int i=0;i<length;i++){out[i]=source[i%source.length];}return out;
	}
	static void write(String path,String[] names,byte[]... seq){
		check(names.length==seq.length,"Fixture names/bases must pair one-to-one");final ByteStreamWriter out=writer(path);
		try{for(int i=0;i<seq.length;i++){MergeRiboSelectionTest.fasta(out,names[i],seq[i]);}}finally{finish(out);}
	}
	static ArrayList<Read> load(String path){return MergeRiboSelectionTest.load(path);}
	static ByteStreamWriter writer(String path){return MergeRiboSelectionTest.writer(path);}
	static void finish(ByteStreamWriter... writers){MergeRiboSelectionTest.finish(writers);}
	static void check(boolean ok,String reason){MergeRiboSelectionTest.check(ok,reason);}
}
