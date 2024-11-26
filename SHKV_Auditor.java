package kv_whole;//package kv_whole;

import com.google.common.collect.Ordering;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.math3.distribution.BetaDistribution;

import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.IntStream;
//import javafx.util.Pair;

//todo: only the output with ue formatted is need to collected

public class SHKV_Auditor {
	int iter = 10000000;
	// 0: AppData    1: JData     2: gaussian    3: powerlaw     4: amazon    5: cloting
	int dataset_index = 4;   // 0-5
	int[] datasets =  {11,    12,    21,    22,    31,      32};
	int dataset = datasets[dataset_index];
	Random random = new Random();
	//    private HashMap<Integer, Double>[] kvs;
	HashMap<Integer, Double> kvs1;
	HashMap<Integer, Double>[] kvs;
	ArrayList<HashMap<ArrayKey, Integer>>[] frequencyMap = new ArrayList[2];
	ArrayList<HashMap<ArrayKey, Integer>>[] frequencyMap2 = new ArrayList[2];
	ArrayList<HashMap<ArrayKey, Integer>>[] small_frequencyMap = new ArrayList[2];
	//	ArrayList<HashMap<ArrayKey, Integer>>[] frequencyMap_with_key = new ArrayList[2];
	ArrayList<HashMap<ArrayKey, Integer>>[] dist_1 = new ArrayList[2];
	ArrayList<HashMap<ArrayKey, Integer>>[] frequency4mean = new ArrayList[2];
	double[] empirical_eps_between = new double[iter];
	double[] empirical_eps_mean = new double[iter];
//	Processor processor;
	double[] perce = { 0.5, 0.5, 0.5};

	int n;
	int d;
	double user_id;


	int thre = 0;

	int[][] boolof1 = new int[2][iter];

	int g_of_total = 8;

	double[] empirical_eps_uehis = new double[iter];


	int bound = 2 ;
	HashMap<Integer, Double> kvs_ori = new HashMap<>();
	int[] binNum = {4};
	int[][][] perturb_st_predict = new int[2][2][binNum[0]];
	int[][] hier = new int[2][binNum[0]];
	int[][] hier2 = new int[2][binNum[0]];
	int[] perturb_encoding_tmp = new int[binNum[0]];

	double[] boundary = new double[binNum[0]];
	double[][] assignedMean = new double[2][binNum[0]  / 2];
	//    double[] boundary = {-1, -0.5, -0.5, 0.0, 0.0, 0.5, 0.5, 1.0};
	int[][] ptb = new int[2][iter];
	List<Map.Entry<Integer, Double>> kvlist = new ArrayList<>(kvs_ori.entrySet());
	double[] mean_whole = new double[]{0.0, 0.0, 0.0, 0.0};
	double[] count_whole = new double[binNum[0] / 2]; //{0.0, 0.0, 0.0, 0.0};
	public SHKV_Auditor() throws Exception{

		this.d = 1;
		this.n = iter;

		histogramTest();

	}


	public void histogramTest(){
		double[] eps = {0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 4, 6, 8, 10, 16};
		DecimalFormat df = new DecimalFormat("0.############");

		HashMap<Integer, Double>[] kvs = new HashMap[2];
		int method = 0;
		if(method==0)
			System.out.println("\n\t-------------  Histogram Test: PrivKVM*-GVPP: -------------");
		else if(method==1)
			System.out.println("\n\t-------------  Histogram Test: PrivKVM*-SW: -------------");
		else if(method==2)
			System.out.println("\n\t-------------  Histogram Test: PrivKVM*-PM: -------------");
		boundary[0] = -1;
		boundary[binNum[0]-1] = 1;
		double step = (double)2 / (binNum[0]/2);
		for(int i = 1; i < binNum[0] - 1; i = i + 2){
			boundary[i] = boundary[i - 1] + step;
			boundary[i + 1] = boundary[i];
		}
		int meanindex = 0;
		for(int i = 0; i < binNum[0] - 1; i = i + 2){
			assignedMean[0][meanindex] = (boundary[i] + boundary[i + 1])/2;
			assignedMean[1][meanindex] = (boundary[i] + boundary[i + 1])/2;
			meanindex++;
		}
		for(int i=0; i<binNum.length; i++){
			System.out.println("\n******************** bin_num = " + binNum[i] + " ********************");
			for(int j=0; j<eps.length; j++){
//				thre = threshold[j];
				double average_mse = 0;
				double[] histogram_per = null;



				int index = 0;

				meanindex = 0;
				for(int imean = 0; imean < binNum[0] - 1; imean = imean + 2){
					assignedMean[0][meanindex] = (boundary[imean] + boundary[imean + 1])/2;
					assignedMean[1][meanindex] = (boundary[imean] + boundary[imean + 1])/2;
					meanindex++;
				}

				double[] histogram_real = new double[binNum[0] / 2];

				histogram_real[binNum[0]/2-1] = iter/2;
				for(int ih = 0; ih < binNum[0] / 2 - 1; ih++){
					histogram_real[ih] = 0.0;
				}

				for(int iii = 0; iii < 2; iii++){
					frequencyMap[iii] = new ArrayList<>();
					frequencyMap2[iii] = new ArrayList<>();


					frequency4mean[iii] = new ArrayList<>();
					frequencyMap[iii].add(new HashMap<>());
					frequencyMap2[iii].add(new HashMap<>());

					frequency4mean[iii].add(new HashMap<>());
				}

				meanindex = 0;
				assignedMean = new double[2][binNum[0] / 2];
				for(int imean = 0; imean < binNum[0] - 1; imean = imean + 2) {
					assignedMean[0][meanindex] = (boundary[imean] + boundary[imean + 1]) / 2;
					assignedMean[1][meanindex] = (boundary[imean] + boundary[imean + 1]) / 2;
					meanindex++;
				}


//						get_histogram_real(binNum[i], index);
				for(int user_id = 0; user_id < 5; user_id++){
					if(method==0) {
						histogram_per = histogram_PrivKVM(user_id, kvlist, eps[j], binNum[i], user_id);
					}
					double mse = 0;
					for(int l=0; l<binNum[i]/2; l++)
						mse += Math.pow(histogram_real[l]-histogram_per[l], 2)/(binNum[i]/2);

					average_mse += mse/d;

				}
			}
		}
	}




	public double[] get_histogram(int binNum, double[] value){
		double[] count = new double[binNum];
		double[] boundary = getBoundary(binNum);
		for(int i=0; i<value.length; i++){
			if(value[i] == Double.MAX_VALUE)
				continue;
			if(value[i]>=boundary[0] && value[i]<=boundary[1])  // in the first interval, the left boundary should be included
				count[0]++;
			for(int j=1; j<boundary.length-1; j++){
				if(value[i]>boundary[j] && value[i]<=boundary[j+1])
					count[j]++;
			}
		}
		return count;
	}

	// set the boundary with binNum bins, in ascending order
	public double[] getBoundary(int binNum){
		double[] boundary = new double[binNum+1];
		boundary[0] = -1;
		boundary[binNum] = 1;
		for(int i=1; i<binNum; i++)
			boundary[i] = boundary[0] + 2.0/binNum*i;
		return boundary;
	}



	public double PiecewiseMechanismPerturbation(double value, double epsilon){
		double res = 0;
		double C = (Math.exp(epsilon/2)+1)/(Math.exp(epsilon/2)-1);
		double p = (Math.exp(epsilon) - Math.exp(epsilon/2)) / (2*Math.exp(epsilon/2) + 2);
		double l = (C+1)/2*value - (C-1)/2;
		double r = l + C - 1;

		Random random = new Random();
		if(random.nextDouble()<p){   // sample a value from [l,r]
			res = getRandomNum(l,r);
		}else{
			res = getRandomNum(-1*C, C);
			while(res<=r && res>=l)
				res = getRandomNum(-1*C, C);
		}
		res /= C;
		return res;
	}

	public static double getRandomNum(double l, double r){
		Random random = new Random();
		double rnd = l + (r-l)*random.nextDouble();
		return rnd;
	}


	public double SquarWaveMechanismPerturbation(double value, double epsilon){
		double res = 0;

		value = (value+1)/2.0;
		double temp = Math.exp(epsilon);
		double b = (epsilon*temp-temp+1)/(2*temp*(temp-1-epsilon));
		double p = temp/(2*b*temp+1);
		double q = 1.0/(2*b*temp+1);

		Random random = new Random();
		if(random.nextDouble()<p){
			res = value - b + (2*b) * random.nextDouble();
		}else{
			res = -b + (2*b+1) * random.nextDouble();
			while(res>value-b && res<value+b)
				res = -b + (2*b+1) * random.nextDouble();
		}

		res = 2*(res+b)/(2*b+1)-1;

		return res;
	}


	public double[] getAllValuesAfterKeyPerturbation(int index, Map.Entry<Integer, Double> kvs, double epsilon, int binNum, int times, int iteration, double ratio_p){
		double p = Math.exp(epsilon)/(1.0+Math.exp(epsilon));
		if(times > 0)
			p = 0.5;

		double[] cnt_v = new double[binNum];

		double v = kvs.getValue();

		double[] value = new double[iteration];
		int cnt = 0; //是根据user的个数进行计数的

		for(int i = 0; i < iteration; i++){
			if(random.nextDouble() < ratio_p){
				boolof1[index][i] = 1;
			}else
				boolof1[index][i] = 0;

			if((boolof1[index][i]== 1)) {
				if (random.nextDouble() < p) {
					value[cnt++] = v;  // 1->1
					ptb[index][i] = 1;
				}
				else {
					value[cnt++] = Double.MAX_VALUE;  // 1->0
					//			}else{
					ptb[index][i] = 0;
				}
			}else{
				if (random.nextDouble() < p) {  // 0->0
					value[cnt++] = Double.MAX_VALUE;
					ptb[index][i] = 1;
				}
				else {
					int idx = random.nextInt(binNum/2);
					value[cnt++] = assignedMean[index][idx];  //哦这块是对它的扰动 就是用mean的部分
					ptb[index][i] = 1;
				}
			}
		}
		return value;

	}//todo：没懂 这个value的赋值是？ 所有的value是只看哪一位？


	public double calculate_cv(ArrayList<HashMap<ArrayKey, Integer>> input_list, int threshold){
		int b_nun = 0;
		double sum_befor = 0.0;
		for (int num : input_list.get(0).values()) {
			if(num < threshold)
				continue;
			sum_befor += num;
			b_nun++;
		}
		sum_befor = sum_befor / b_nun;
		// Step 2: 计算方差
		double variance_sort = 0.0;
		for (int num : input_list.get(0).values()) {
			if(num < threshold)
				continue;
			variance_sort += Math.pow(num - sum_befor, 2);
		}
		variance_sort /= b_nun;
		variance_sort = Math.sqrt(variance_sort);
		double cov_sort = variance_sort / sum_befor;
		return cov_sort;
	}



	public double calculate_cv_ori(ArrayList<HashMap<ArrayKey, Integer>> input_list, int threshold){
		int b_nun = 0;
		double sum_befor = 0.0;
		for (int num : input_list.get(0).values()) {
			if(num < threshold)
				continue;
			sum_befor += num;
			b_nun++;
		}
		sum_befor = sum_befor / b_nun;
		double variance_sort = 0.0;
		for (int num : input_list.get(0).values()) {
			if(num < 0)
				continue;
			variance_sort += Math.pow(num - sum_befor, 2);
		}
		variance_sort /= b_nun;

		variance_sort = Math.sqrt(variance_sort);
		double cov_sort = variance_sort / sum_befor;
		return cov_sort;
	}

	public double[] calculate_threshold(ArrayList<HashMap<ArrayKey, Integer>> change, double ori_cv, double eps, String str1){
		double cal_cv = 0.0;
		int num_thresh = 5;
		if(eps <= 0.2)
			cal_cv = 3;

		else if(eps <= 0.4)
			cal_cv = 3;
		else
			cal_cv = 6;

		double sum2 = 0.0;
		double variance2 = 0.0;
		int threshold = 0;
		TreeMultimap<Integer, String> s1 = TreeMultimap.create(
				Ordering.natural().reverse(), // 降序排列键
				Ordering.natural()
		);

		for (HashMap<ArrayKey, Integer> entry : change) {
			for(ArrayKey key : entry.keySet()) {
				s1.put(change.get(0).get(key), Arrays.toString(key.array));
			}
		}
		List<ArrayKey> intersection = new ArrayList<>();

//		intersection = new ArrayList<>();
		int num_tree = 0;
		double cov_tree1 = 0.0;
		int num_1 = 0;
//		for(int j = 0; j < s1.size(); j++) {
		for(Map.Entry<Integer, String> entry : s1.entries()) {
			int num = entry.getKey();
			if(num == 1){
				num_1++;
			}
			if(num < num_thresh){
				break;
			}
			num_tree++;
			sum2 += num;
			sum2 /= num_tree;
			int t_tmp = 0;
			variance2 = 0.0;
			for (Map.Entry<Integer, String> entry_t : s1.entries()) {
				if (t_tmp < num_tree) {
					int num_t = entry_t.getKey();
					variance2 += Math.pow(num_t - sum2, 2);
					t_tmp++;
				} else {
					break;
				}
			}
			variance2 /= num_tree;
			variance2 = Math.sqrt(variance2);
			cov_tree1 = variance2 / sum2;
			if (cov_tree1 >= cal_cv) {
				num_tree--;
//				System.out.println("num : " + num);
				break;
			} else {
				sum2 *= num_tree;
				String str = new String(entry.getValue());
				str = str.substring(1, str.length() - 1);
				String[] stringArray = str.split(",\\s*");
				int[] intArray = new int[stringArray.length];
				for (int i = 0; i < stringArray.length; i++) {
					intArray[i] = Integer.parseInt(stringArray[i]);
				}
				ArrayKey tmp_a = new ArrayKey(intArray);
				if(!intersection.contains(tmp_a)) {
					intersection.add(tmp_a);
				}
			}
		}
//		}

		int numtmp = 1;
		int flag = 0;
		for(Map.Entry<Integer, String> entry : s1.entries()) {
			threshold = entry.getKey();
			if(numtmp == num_tree) {
				break;
			}
			numtmp++;
			flag = 1;
		}

		numtmp = 0;
		double sum = 0.0;
		for (Map.Entry<Integer, String> entry : s1.entries()) {
			if(numtmp == num_tree)
				break;
			double value = entry.getKey();
			sum += value;
		}
		double mean = sum / num_tree;
		numtmp = 0;
		for (Map.Entry<Integer, String> entry : s1.entries()) {
			if(numtmp == num_tree)
				break;
			double value = entry.getKey();
			sum += Math.pow(value - mean, 2);
		}


		double[] result = new double[2];
		result[0] = threshold;
		result[1] = 0;

		return result;
	}

	public double[] histogram_PrivKVM(int index,  List<Map.Entry<Integer, Double>> kvs, double epsilon, int binNum, int times){
		double[] count_per = null;
		DecimalFormat df = new DecimalFormat("0.############");
		for(int i = 0; i < 2; i++){
			frequencyMap[i] = new ArrayList<>();
			frequencyMap2[i] = new ArrayList<>();
			small_frequencyMap[i] = new ArrayList<>();

			frequencyMap[i].add(new HashMap<>());
			frequencyMap2[i].add(new HashMap<>());
			small_frequencyMap[i].add(new HashMap<>());

		}

		String uebeforformatted = new String();
		hier = new int[2][binNum];
		hier2 = new int[2][binNum];

		System.out.println("assign : " + assignedMean[0][0]);
		kvs_ori.put(2,  assignedMean[0][0]);
		kvs_ori.put(1,  assignedMean[0][0]);
		kvs.addAll(new HashMap<>(kvs_ori).entrySet());
		for(int i = 0; i < binNum / 2; i++){
			assignedMean[1][i] = assignedMean[0][i];
		}
		int perce_ratio = 1;
		perturb_st_predict = new int[2][2][binNum];
		count_per = histogram_UE(1, binNum, epsilon, kvs, times, perce_ratio);

		int key_k = kvs.get(0).getKey();
		kvs.clear();
		kvs_ori.put(2, -1.0);
		kvs_ori.put(1, -0.5);
		kvs.addAll(new HashMap<>(kvs_ori).entrySet());


		count_per = histogram_UE(0, binNum, epsilon, kvs, times, perce_ratio);
		kvs.clear();

		System.out.println("eps" + epsilon);
		empirical_eps_mean = new double[iter];
		double[] eps_hier = new double[binNum];

		int test_statistic = 0;
		System.out.println(perturb_st_predict[0][0][test_statistic] );
		double[] bound1_pre = clopperPearson(perturb_st_predict[0][0][test_statistic] , iter, 0.05);
		double[] bound2_pre = clopperPearson(perturb_st_predict[1][0][test_statistic], iter, 0.05);
		double log1_pre = Math.abs(Math.log(bound1_pre[0] / bound2_pre[1]));
		double log2_pre = Math.abs(Math.log(bound1_pre[1] / bound2_pre[0]));
		System.out.println("ac的 ： " +  Math.max(log2_pre, log1_pre));

		System.out.println(perturb_st_predict[0][1][test_statistic] );
		bound1_pre = clopperPearson(perturb_st_predict[0][1][test_statistic] , iter, 0.05);
		bound2_pre = clopperPearson(perturb_st_predict[1][1][test_statistic], iter, 0.05);
		log1_pre = Math.abs(Math.log(bound1_pre[0] / bound2_pre[1]));
		log2_pre = Math.abs(Math.log(bound1_pre[1] / bound2_pre[0]));
		System.out.println("ac的 ： " +  Math.max(log2_pre, log1_pre));

		for(int k = 0; k < binNum; k++){
			double[] bound1 = clopperPearson(hier[0][k], iter/2, 0.05);
			double[] bound2 = clopperPearson(hier[1][k], iter/2, 0.05);
			double eps1 = Math.abs(Math.log(bound1[0] / bound2[1]));
			double eps2 = Math.abs(Math.log(bound1[1] / bound2[0]));

			eps_hier[k] = Math.max(eps2, eps1);

		}
		double max_eps_hier = Arrays.stream(eps_hier).max().getAsDouble();
		System.out.println("hier : " +max_eps_hier );

		eps_hier = new double[binNum];
		for(int k = 0; k < binNum; k++){
			double[] bound1 = clopperPearson(hier2[0][k], iter/2, 0.05);
			double[] bound2 = clopperPearson(hier2[1][k], iter/2, 0.05);
			double eps1 = Math.abs(Math.log(bound1[0] / bound2[1]));
			double eps2 = Math.abs(Math.log(bound1[1] / bound2[0]));
			eps_hier[k] = Math.max(eps2, eps1);

		}
		max_eps_hier = Arrays.stream(eps_hier).max().getAsDouble();
		System.out.println("hier2 : " +max_eps_hier );

		for(int k = 0; k < binNum; k++){
			double[] bound1 = clopperPearson(hier[0][k], iter/2, 0.05);
			double[] bound2 = clopperPearson(hier2[0][k], iter/2, 0.05);
			double eps1 = Math.abs(Math.log(bound1[0] / bound2[1]));
			double eps2 = Math.abs(Math.log(bound1[1] / bound2[0]));
			eps_hier[k] = Math.max(eps2, eps1);
		}
		max_eps_hier = Arrays.stream(eps_hier).max().getAsDouble();
		System.out.println("hier normal : " +max_eps_hier );

//		List<ArrayKey> intersection_small = new ArrayList<>();
		HashMap<ArrayKey, Integer> tmp = new HashMap<>();

		HashMap<ArrayKey, Integer> tmp2 = new HashMap<>();
		ArrayList<HashMap<ArrayKey, Integer>>[] inter_frequencyMap = new ArrayList[2];

		List<ArrayKey> calculate = new ArrayList<>();
		ArrayList<ArrayKey> intersection = new ArrayList<>();



		int id = 0;

		int min_add = 0;


		int thre1 = 0, thre2 = 0;
		tmp = frequencyMap2[0].get(0);

		inter_frequencyMap[0] = new ArrayList<>();
		inter_frequencyMap[1] = new ArrayList<>();
		inter_frequencyMap[0].add(new HashMap<>());
		inter_frequencyMap[1].add(new HashMap<>());
		intersection = new ArrayList<>();
		//todo 直接扰动完的inter
		for (ArrayKey intmp : tmp.keySet()) {
			tmp2 = frequencyMap[1].get(0);
			if (tmp2.containsKey(intmp)) {
				intersection.add(intmp);
				inter_frequencyMap[0].get(0).put(intmp, tmp.get(intmp));
				inter_frequencyMap[1].get(0).put(intmp, tmp2.get(intmp));
				//							}
			}
		}
		double cv_cal1 = calculate_cv_ori(inter_frequencyMap[0], 0);
		double cv_cal2 = calculate_cv_ori(inter_frequencyMap[1], 0);
//		System.out.println("ue formatted 原来的cv： " + cv_cal1 + " " + cv_cal2 );

		int[] tmp_his = new int[iter];
		int ind = 0;
		for(ArrayKey thm: frequencyMap[0].get(0).keySet()){
			int inte = frequencyMap[0].get(0).get(thm);
			tmp_his[ind] = inte;
			ind++;
		}



		double[] ans1 = calculate_threshold(inter_frequencyMap[0], cv_cal1, epsilon, "normal");
		thre1 = (int)ans1[0];
//		thre1 = 1;
		double[] ans2 = calculate_threshold(inter_frequencyMap[1], cv_cal2, epsilon, "normal");
		thre2 = (int)ans2[0];
//		thre2 = 1;
		double cv_1 = calculate_cv(inter_frequencyMap[0], thre1);
		double cv_2 = calculate_cv(inter_frequencyMap[1], thre2);
		double log1 = 0, log2 = 0;
		empirical_eps_uehis = new double[iter];
		calculate = new ArrayList<>();

		for (ArrayKey key : intersection) {

			double[] bound_before1 = new double[2];
			double[] mean1 = new double[2];
			double[] mean2 = new double[2];

			double[] bound1 = new double[2];
			double[] bound2 = new double[2];
			if (frequencyMap2[1].get(0).get(key)  < thre2|| frequencyMap[1].get(0).get(key) < thre1) {
				continue;
			} else {
				bound1 = clopperPearson((int) frequencyMap2[1].get(0).get(key), iter, 0.05);
				bound2 = clopperPearson((int) frequencyMap[1].get(0).get(key), iter, 0.05);
			}
			log1 = Math.abs(Math.log(bound1[0] / bound2[1]));
			log2 = Math.abs(Math.log(bound1[1] / bound2[0]));
			empirical_eps_uehis[id] = Math.max(log1, log2);
			empirical_eps_uehis[id] = Math.abs(empirical_eps_uehis[id]);
			calculate.add(key);
			id = id + 1;

		}

		id = 0;

		uebeforformatted = df.format(Arrays.stream(empirical_eps_uehis).max().getAsDouble());
		OptionalInt maxIndex = IntStream.range(0, empirical_eps_uehis.length)
				.reduce((i, j) -> empirical_eps_uehis[i] > empirical_eps_uehis[j] ? i : j);

		int index_max = maxIndex.getAsInt();
		int tmp_index = 0;
		ArrayKey maxkey = new ArrayKey(perturb_encoding_tmp);
		for(ArrayKey key : calculate){
			maxkey = key;

			if(tmp_index== index_max){
				break;
			}
			tmp_index++;
		}
		double[] empirical_eps_uehis2 =new double[iter];
		for (ArrayKey key : intersection) {

			double[] bound_before1 = new double[2];
			double[] mean1 = new double[2];
			double[] mean2 = new double[2];

			double[] bound1 = new double[2];
			double[] bound2 = new double[2];
			if (frequencyMap2[0].get(0).get(key)  < thre2|| frequencyMap[0].get(0).get(key) < thre1) {
				continue;
//				empirical_eps_uehis[id] = 0.0;
			} else {
				bound1 = clopperPearson((int) frequencyMap2[0].get(0).get(key), iter, 0.05);
				bound2 = clopperPearson((int) frequencyMap[0].get(0).get(key), iter, 0.05);
			}
			log1 = Math.abs(Math.log(bound1[0] / bound2[1]));
			log2 = Math.abs(Math.log(bound1[1] / bound2[0]));
			empirical_eps_uehis2[id] = Math.max(log1, log2);
			empirical_eps_uehis2[id] = Math.abs(empirical_eps_uehis2[id]);

			calculate.add(key);
			id = id + 1;

		}
		String uebeforformatted2 = new String();
		uebeforformatted2 = df.format(Arrays.stream(empirical_eps_uehis2).max().getAsDouble());
		id = 0;
		System.out.println("ue formatted " + uebeforformatted +  uebeforformatted2 + " " +thre1 + " " + thre2+ " " + cv_1 + " " + cv_2 + " " +frequencyMap[1].get(0).get(maxkey) + " " + frequencyMap[0].get(0).get(maxkey));

		intersection = new ArrayList<>();
		inter_frequencyMap[0] = new ArrayList<>();
		inter_frequencyMap[1] = new ArrayList<>();
		inter_frequencyMap[0].add(new HashMap<>());
		inter_frequencyMap[1].add(new HashMap<>());

		for (ArrayKey intmp : tmp.keySet()) {
			tmp2 = frequencyMap[0].get(0);
			//						for (int listi2 = 0; listi2 < frequencyMap[1].size(); listi2++) {
			//							Map<ArrayKey, Integer> tmp2 = frequencyMap[1].get(listi2);
			if (tmp2.containsKey(intmp)) {
				intersection.add(intmp);
				inter_frequencyMap[0].get(0).put(intmp, tmp.get(intmp));
				inter_frequencyMap[1].get(0).put(intmp, tmp2.get(intmp));
			}
		}

		id = 0;
		return count_per;
	}


	public double[] histogram_UE(int index, int binNum, double epsilon,   List<Map.Entry<Integer, Double>> kvs, int times, int perce_idx){
		double[] values = new double[iter];
		double[] values2 = new double[iter];
		int iteration = iter;
		if(kvs.size() == 2)
			iteration = iter / 2;
		double perce_ratio =0.0;
		perce_ratio = perce[times % 3];

		if(kvs.size() == 2) {
			values = getAllValuesAfterKeyPerturbation(index, kvs.get(0), epsilon / 2.0, binNum, times, iteration, perce_ratio); //todo 是不是先把数据给扰动了啊？
//        double p = Math.exp(epsilon)/(1.0+Math.exp(epsilon));
			if(times % 3 + 1 < 3){
				perce_ratio = perce[times % 3 + 1];
			}else{
				perce_ratio = perce[times % 3 + 1 - 3];
			}
//			perce_ratio = perce[2 - times % 3];
			values2 = getAllValuesAfterKeyPerturbation(index, kvs.get(1), epsilon / 2.0, binNum, times, iteration, perce_ratio); //todo 是不是先把数据给扰动了啊？
		}
		else
			values = getAllValuesAfterKeyPerturbation(index, kvs.get(0), epsilon / 2.0, binNum, times, iteration, perce_ratio); //todo 是不是先把数据给扰动了啊？

		double[] cnt_v = new double[binNum];
		double[] c_m = new double[binNum];

		//在这里就可以加入对values分布的统计了吧 不行 还得看那个扰动之后的unary的情况 但是key的话要考虑加入进来
		double[] noisyCount = new double[binNum];
//        double[] noisy_Count_fuben = new double[binNum];
		double[] noisyCount_with_key = new double[binNum + 1];
		int[] perturb_encoding_with_key = new int[binNum + 1];

		for(int i=0; i<values.length; i++) {
			int[] perturb_encoding = valuePert_UnaryEncoding_histogram_4mean(binNum, epsilon / 2.0/5, values[i], i);
			for (int k = 0; k < noisyCount.length; k++)
				noisyCount[k] += perturb_encoding[k];  //这个就是没有扰校正的那个OUE的情况
			for(int k = 0; k < hier[index].length; k++){
				hier[index][k] += perturb_encoding[k]; //这个就是数据的
			}
			ArrayKey key = new ArrayKey(perturb_encoding);
			if(index ==0){
				if (!frequencyMap2[0].get(0).containsKey(key)) {
					frequencyMap2[0].get(0).put(key, 1);
				} else {
					int value = frequencyMap2[0].get(0).get(key);
					frequencyMap2[0].get(0).put(key, value + 1);
				}
			}else if(index ==1) {
//			ArrayKey key = new ArrayKey(mix_perturb_encoding);
				if (!frequencyMap[0].get(0).containsKey(key)) {
					frequencyMap[0].get(0).put(key, 1);
				} else {
					int value = frequencyMap[0].get(0).get(key);
					frequencyMap[0].get(0).put(key, value + 1);
				}
			}
			ArrayList<Integer> perturb_st_predict_idx = new ArrayList<>();
			for(int k = 0; k < noisyCount.length; k++){
				if(perturb_encoding[k] == 1)
					perturb_st_predict_idx.add(k);
			}

			int predict_idx = 0;
			if(perturb_st_predict_idx.size() == 0){
				predict_idx = random.nextInt(noisyCount.length);
			}else {
				predict_idx = perturb_st_predict_idx.get(random.nextInt(perturb_st_predict_idx.size()));
			}
			perturb_st_predict[index][0][predict_idx] += 1;
		}
		if(kvs.size() == 2) {
			for (int i = 0; i < values2.length; i++) {
				int[] perturb_encoding = valuePert_UnaryEncoding_histogram_4mean(binNum, epsilon / 2.0/5, values2[i], i);
				for (int k = 0; k < noisyCount.length; k++)
					noisyCount[k] += perturb_encoding[k];
				for (int k = 0; k < hier2[index].length; k++) {
					hier2[index][k] += perturb_encoding[k]; //这个就是数据的 第二个扰动部分的
				}
				ArrayList<Integer> perturb_st_predict_idx = new ArrayList<>();
				for(int k = 0; k < noisyCount.length; k++){
					if(perturb_encoding[k] == 1)
						perturb_st_predict_idx.add(k);
				}

				ArrayKey key = new ArrayKey(perturb_encoding);
//			ArrayKey key = new ArrayKey(mix_perturb_encoding);

				if(index ==0){
					if (!frequencyMap2[1].get(0).containsKey(key)) {
						frequencyMap2[1].get(0).put(key, 1);
					} else {
						int value = frequencyMap2[1].get(0).get(key);
						frequencyMap2[1].get(0).put(key, value + 1);
					}
				}else if(index ==1) {
//			ArrayKey key = new ArrayKey(mix_perturb_encoding);
					if (!frequencyMap[1].get(0).containsKey(key)) {
						frequencyMap[1].get(0).put(key, 1);
					} else {
						int value = frequencyMap[1].get(0).get(key);
						frequencyMap[1].get(0).put(key, value + 1);
					}
				}
				int predict_idx = 0;
				if(perturb_st_predict_idx.size() == 0){
					predict_idx = random.nextInt(noisyCount.length);
				}else {
					predict_idx = perturb_st_predict_idx.get(random.nextInt(perturb_st_predict_idx.size()));
				}
				perturb_st_predict[index][1][predict_idx] += 1;
			}

		}
		int numberOfValues = 0;
		for(double v : values)
			if(v<Double.MAX_VALUE)
				numberOfValues++; //统计每一个value它的个数
		if(index == 0)
			for(double v : values2)
				if(v<Double.MAX_VALUE)
					numberOfValues++;
		calibrate(noisyCount, 2, numberOfValues, epsilon/2.0/5);  //这个是单独校正count todo：和calibrate_oue相同的  // epsilon for Unary Encoding with sensitivity of 2
		//这个也是校正的noisy count啊


		int j = 0;
		for(int i = 0; i < noisyCount.length ; i = i+ 2){
			assignedMean[index][j] = (boundary[i+1] * noisyCount[i+1] + boundary[i]*noisyCount[i]) / (noisyCount[i] + noisyCount[i + 1]) ;
			j = j + 1;
		}
		//todo：确实是聚合部分的校正 是每个bucket的校正
		int counti = 0;
		for(int i=0; i<noisyCount.length; i = i + 2) { //第二次校正 应该是聚合部分的校正吧
//			noisyCount[i] = rangeCountCalibration(noisyCount[i] +noisyCount[i+1] , 0, epsilon, binNum+1, numberOfValues);  // key is perturbed with epsilon
			count_whole[counti] = rangeCountCalibration(noisyCount[i] + noisyCount[i + 1], 0, epsilon / 2.0/5, binNum / 2 + 1, numberOfValues);
			counti = counti + 1;
//			System.out.println("count: " + count_whole[counti - 1]);
		}
		return count_whole;  //这个是用公式17进行校正的
	}


	//todo 但是你这个计算你也只是一个数据去计算啊 你直方图去计算的话只有一个左边界的就算是 去进行扰动的
	// due to key sampling
	public double rangeCountCalibration(double v, int key, double epsilon, int g, int total){
//		double frequency = this.processor.frequency[key];
		double p = Math.exp(epsilon)/(1+Math.exp(epsilon));
//		double res = v/p - this.n*(1-frequency)*(1-p)/(g-1)/p;  //公式17
//		double frequency = (p - 1 + total) / (2 * p - 1);
		double frequency = (double)total / iter;
		frequency = (p - 1 + frequency) / (2 * p - 1);
		double res = iter*v/p/total - iter*(1-frequency)*(1-p)/(g-1)/p;
		if(res<0)
			res = 0;
		if(res>total)
			res = total;
		return res;
	}


	public int valuePert_GRR_histogram(int binNum, double ep, double value){
		if(value == Double.MAX_VALUE)
			return Integer.MAX_VALUE;
		int idx_bin = rangeDiscretization_histogram(value, binNum);  //看这个值是落到了哪个bucket中的
		return GRR(idx_bin, binNum, ep);  //把这个idx给它进行扰动的
	}


	public static int GRR(int index, int num, double epsilon){
		double exp = Math.exp(epsilon);
		double p = exp/(exp+num-1);
		Random random = new Random();
		if(random.nextDouble() < p)
			return index;
		int idx = random.nextInt(num);
		while(idx==index)
			idx = random.nextInt(num); //哦还可以这样找不同的啊？
		return idx;
	}


	public int rangeDiscretization_histogram(double value, int binNum){
		double[] boundary = getBoundary(binNum);
		if(value>=boundary[0] && value<=boundary[1])
			return 0;
		int index = 1;
		for(int i=1; i<binNum; i++){
			if(value>boundary[i] && value<=boundary[i+1]){
				index = i;
				break;
			}
		}
		return index;
	}

	public int rangeDiscretization_histogram_4mean(double value, int binNum, int in){

		double p = 0.0;
		Random random = new Random();
		if(value>=boundary[0] && value<=boundary[1]){
			p = (boundary[1] - value) / (boundary[1] - boundary[0]);
			if(random.nextDouble() < p) {
//				mean_index[in] = 0;
				return 0;

			}else{
//				mean_index[in] = 1;
				return 1;
			}
		}
		int index = 1;
		for(int i=2; i<binNum; i = i + 2){
			if(value>boundary[i] && value<=boundary[i+1]){
				p = (boundary[i + 1] - value) / (boundary[i+1] - boundary[i]);
				if(random.nextDouble() < p){
//                    index = 2 * i;
					index = i;
				}else{
//                    index = 2 * i + 1;
					index = i+1;
				}
				break;
			}
		}
//		mean_index[in] = index;
		return index;
	}


	public static double[] valueCountCalibration(int[] perturb_index, int num, double ep){
		double[] count = new double[num];
		int cnt = 0;
		for(int i=0; i<perturb_index.length; i++){
			if(perturb_index[i]<Integer.MAX_VALUE){
				count[perturb_index[i]]++; //每个perturb的index的个数
				cnt++;  //总的样本个数
			}
		}
		calibrate(count, num, cnt, ep);
		return count;
	}

	public static void calibrate(double[] count, int k, int N, double ep){
		double exp = Math.exp(ep);
		double p = exp / (k-1+exp);
		for(int i=0; i<count.length; i++){
			count[i] = (int)Math.round( N*(p-1)/(k*p-1.0) + count[i]*(k-1.0)/(k*p-1.0) );  //todo：这个公式真的不知道是哪里来的啊
//			count[i] = (2 *(exp + 1)* count[i] - 2 * N ) / (exp - 1);
			if(count[i]>N)
				count[i] = N;
			if(count[i]<0)
				count[i] = 0;
		}
	}


	public int[] valuePert_UnaryEncoding_histogram_4mean(int binNum, double ep, double value, int in){
		if(value == Double.MAX_VALUE)
			return new int[binNum];
		int idx = rangeDiscretization_histogram_4mean(value, binNum, in); // 1->left, 2->right, 0->others  //让他转换成直方图的形式
		return UnaryEncoding_4mean(idx, binNum , ep);  //返回的是编码
	}

	public int[] valuePert_UnaryEncoding_histogram(int binNum, double ep, double value){
		if(value == Double.MAX_VALUE)
			return new int[binNum];
		int idx = rangeDiscretization_histogram(value, binNum); // 1->left, 2->right, 0->others  //让他转换成直方图的形式
		return UnaryEncoding(idx, binNum, ep);  //返回的是编码
	}
	public static double[] clopperPearson(int k, int n, double alpha) {
		double lowerBound, upperBound;


		BetaDistribution betaLower = new BetaDistribution(k, n - k + 1);
		lowerBound = betaLower.inverseCumulativeProbability(alpha / 2);


		BetaDistribution betaUpper = new BetaDistribution(k + 1, n - k);
		upperBound = betaUpper.inverseCumulativeProbability(1 - alpha / 2);


		return new double[]{lowerBound, upperBound};
	}

	public static int[] UnaryEncoding_4mean(int index, int num, double ep){
		double exp = Math.exp(ep);
		double p = exp/(exp+1);
		int[] encoding = new int[num];
		Random random = new Random();
//        int index_4mean = index / 2;
		for(int i=0; i<num; i++){
			if(i==index){
				if(random.nextDouble()<=p) //生成一个随机值 (0,1]范围之内的
					encoding[i] = 1;
			}else{
				if(random.nextDouble()>p)
					encoding[i] = 1;
			}
		}
		return encoding;
	}

	public static int[] UnaryEncoding(int index, int num, double ep){
		double exp = Math.exp(ep);
		double p = exp/(exp+1);
		int[] encoding = new int[num];
		Random random = new Random();
		for(int i=0; i<num; i++){
			if(i==index){
				if(random.nextDouble()<=p) //生成一个随机值 (0,1]范围之内的
					encoding[i] = 1;
			}else{
				if(random.nextDouble()>p)
					encoding[i] = 1;
			}
		}
		return encoding;
	}

	public static void main(String[] args) throws Exception{
		new SHKV_Auditor();
	}
}

