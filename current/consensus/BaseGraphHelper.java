package consensus;

public class BaseGraphHelper {

	public static void initForScoring(BaseGraph bg){
		bg.calcProbs();
	}

	public static void setCounts(BaseGraph bg, int pos, int a, int c, int g, int t, int del){
		BaseNode rnode=bg.ref[pos];
		rnode.acgtCount[0]=a;
		rnode.acgtCount[1]=c;
		rnode.acgtCount[2]=g;
		rnode.acgtCount[3]=t;
		rnode.countSum=a+c+g+t;
		rnode.weightSum=a+c+g+t;
		rnode.acgtWeight[0]=a;
		rnode.acgtWeight[1]=c;
		rnode.acgtWeight[2]=g;
		rnode.acgtWeight[3]=t;
		bg.del[pos].countSum=del;
		bg.del[pos].weightSum=del;
	}

	public static int[] getCounts(BaseGraph bg, int pos){
		BaseNode rnode=bg.ref[pos];
		int del=bg.del[pos].countSum;
		return new int[]{rnode.acgtCount[0], rnode.acgtCount[1],
			rnode.acgtCount[2], rnode.acgtCount[3], del};
	}

	public static int length(BaseGraph bg){
		return bg.ref.length;
	}
}
