package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.io.BufferedReader;
import java.util.TreeMap;
import java.util.HashSet;

/** Counts-only comparison of existing grid criteria with S/N; never selects ranks.
 * Usage: java -ea prot.MarkerYieldSummary RAW_TSV GRID_TSV
 * Author/reviewer: Yoimiya. */
public class MarkerYieldSummary {
    static final int[] CUTS={95,90,85,80,75,70};
    static class Group {
        String kind; int n; int strict=-1, floor=-1;
        int[] sn=new int[CUTS.length];
        HashSet<Integer> ranks=new HashSet<Integer>();
        Group(String k,int size){kind=k;n=size;}
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=2){throw new IllegalArgumentException("RAW_TSV GRID_TSV required");}
        TreeMap<String,Group> groups=new TreeMap<String,Group>();
        try(BufferedReader in=Files.newBufferedReader(Paths.get(args[1]),StandardCharsets.UTF_8)){
            for(String line;(line=in.readLine())!=null;){
                if(line.startsWith("#")||line.isEmpty()){continue;}
                String[] f=line.split("\t",-1);
                if(f.length!=6){throw new IllegalArgumentException("Grid width");}
                int n=Integer.parseInt(f[2]), count=Integer.parseInt(f[5]);
                if(n<1||count<0){throw new IllegalArgumentException("Negative grid count");}
                Group g=groups.get(f[0]);
                if(g==null){g=new Group(f[1],n);groups.put(f[0],g);}
                if(g.n!=n||!g.kind.equals(f[1])){throw new IllegalArgumentException("Inconsistent group");}
                double p=Double.parseDouble(f[3]),s=Double.parseDouble(f[4]);
                if(p==0.95&&s==0.90){if(g.strict!=-1){throw new IllegalArgumentException("Duplicate strict cell");}g.strict=count;}
                if(p==0.70&&s==0.70){if(g.floor!=-1){throw new IllegalArgumentException("Duplicate floor cell");}g.floor=count;}
            }
        }
        try(BufferedReader in=Files.newBufferedReader(Paths.get(args[0]),StandardCharsets.UTF_8)){
            for(String line;(line=in.readLine())!=null;){
                if(line.startsWith("#")||line.isEmpty()){continue;}
                String[] f=line.split("\t",-1);
                if(f.length!=6){throw new IllegalArgumentException("Raw width");}
                Group g=groups.get(f[0]);
                int n=Integer.parseInt(f[2]),r=Integer.parseInt(f[3]);
                int p=Integer.parseInt(f[4]),s=Integer.parseInt(f[5]);
                if(g==null||n!=g.n||!f[1].equals(g.kind)||s<0||p<s||p>n||!g.ranks.add(r)){
                    throw new IllegalArgumentException("Invalid or duplicate raw row: "+line);
                }
                for(int i=0;i<CUTS.length;i++){if(100L*s>=(long)CUTS[i]*n){g.sn[i]++;}}
            }
        }
        System.out.print("#group\tkind\tn_orgs\tP95_SgivenP90\tP70_SgivenP70");
        for(int cut:CUTS){System.out.print("\tSofN"+cut);}
        System.out.println();
        for(java.util.Map.Entry<String,Group> entry:groups.entrySet()){
            Group g=entry.getValue();
            if(g.strict<0||g.floor<0){throw new IllegalArgumentException("Missing comparison cells");}
            System.out.print(entry.getKey()+"\t"+g.kind+"\t"+g.n+"\t"+g.strict+"\t"+g.floor);
            for(int count:g.sn){System.out.print("\t"+count);}
            System.out.println();
        }
    }
}
