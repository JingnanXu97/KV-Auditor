package th_w;
//import ArrayKey;
import org.apache.commons.math3.distribution.BetaDistribution;

import java.util.*;
import java.util.stream.Collectors;

public class the {


    int iter = 100000000;
    int v1 = 0;
    int v2 = 1;
    int binNum = 4;
    double[] epsilon = {0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 4, 6};

    ArrayList<HashMap<List<Double>, Integer>>[] frequencyMap = new ArrayList[2];
    ArrayList<HashMap<List<Double>, Integer>>[] frequencyMap_2bit = new ArrayList[2];
    public the() throws Exception{
       audit();
    }

    public void audit(){
        double[] thresh = {0.5124928545331653,0.5249429746156357,0.5495498751789415,0.596578396322232,0.6769874744974638,0.7823500118298539,0.8156771238097529 ,0.8691255209837462};
        for(int eps = 0; eps < epsilon.length; eps++) {
            double[] pert1 = new double[binNum];
            double[] pert2 = new double[binNum];
            double[] return1 = new double[binNum];
            double[] return2 = new double[binNum];
            double[] emp_eps = new double[binNum];
            frequencyMap[0] = new ArrayList<>();
            frequencyMap[0].add(new HashMap<>());
            frequencyMap[1] = new ArrayList<>();
            frequencyMap[1].add(new HashMap<>());
            frequencyMap_2bit[0] = new ArrayList<>();
            frequencyMap_2bit[0].add(new HashMap<>());
            frequencyMap_2bit[1] = new ArrayList<>();
            frequencyMap_2bit[1].add(new HashMap<>());
            for (int i = 0; i < iter; i++) {
                return1 = THE(v1, -1, 1, epsilon[eps], thresh[eps]);
                List<Double> listKey = Arrays.stream(return1)
                        .boxed()
                        .collect(Collectors.toList());
                List<Double> sele = Arrays.asList(listKey.get(0), listKey.get(1));

                if (!frequencyMap_2bit[0].get(0).containsKey(sele)) {
                    frequencyMap_2bit[0].get(0).put(sele, 1);
                } else {
                    int value = frequencyMap_2bit[0].get(0).get(sele);
                    frequencyMap_2bit[0].get(0).put(sele, value + 1);
                }
                if (!frequencyMap[0].get(0).containsKey(listKey)) {
                    frequencyMap[0].get(0).put(listKey, 1);
                }else{
                    int value = frequencyMap[0].get(0).get(listKey);
                    frequencyMap[0].get(0).put(listKey, value + 1);
                }
                for(int j = 0; j < binNum; j++)
                    pert1[j] += return1[j];
                return2 = THE(v2, -1, 1, epsilon[eps], thresh[eps]);
                listKey = Arrays.stream(return2)
                        .boxed()
                        .collect(Collectors.toList());
                List<Double> sele1 = Arrays.asList(listKey.get(0), listKey.get(1));
//                ArrayKey select = new ArrayKey(sele);
                if (!frequencyMap_2bit[1].get(0).containsKey(sele1)) {
                    frequencyMap_2bit[1].get(0).put(sele1, 1);
                } else {
                    int value = frequencyMap_2bit[1].get(0).get(sele1);
                    frequencyMap_2bit[1].get(0).put(sele1, value + 1);
                }
                if (!frequencyMap[1].get(0).containsKey(listKey)) {
                    frequencyMap[1].get(0).put(listKey, 1);
                }else{
                    int value = frequencyMap[1].get(0).get(listKey);
                    frequencyMap[1].get(0).put(listKey, value + 1);
                }
                for(int j = 0; j < binNum; j++)
                    pert2[j] += return2[j];
            }

            Map<List<Double>, Integer> intersection = new HashMap<>();
            HashMap<List<Double>, Integer> map1 = frequencyMap[0].get(0);
            HashMap<List<Double>, Integer> map2 = frequencyMap[1].get(0);

            for (Map.Entry<List<Double>, Integer> entry : map1.entrySet()) {
                List<Double> key = entry.getKey();
                if (map2.containsKey(key)) {
                    intersection.put(key, Math.min(entry.getValue(), map2.get(key)));
                }
            }


            double[] emp_int = new double[iter];
            int id= 0;
            for (Map.Entry<List<Double>, Integer> entry : intersection.entrySet()) {

                List<Double> listKey = entry.getKey();
                if(frequencyMap[0].get(0).get(listKey) < 1 || frequencyMap[1].get(0).get(listKey) < 1) {
                    continue;
                }
                double[] cp1 = clopperPearson((int)frequencyMap[0].get(0).get(listKey), iter, 0.05);
                double[] cp2 = clopperPearson((int)frequencyMap[1].get(0).get(listKey), iter, 0.05);
                double log_m1 = Math.abs(Math.log(cp1[0] / cp2[1]));
                double log_m2 = Math.abs(Math.log(cp1[1] / cp2[0]));
                emp_int[id] = Math.abs(Math.max(log_m1, log_m2));
                if(Double.isInfinite(emp_int[id])) {
                    emp_int[id] = 0;
                }
                id = id +1;
            }
            System.out.println("hori: " + Arrays.stream(emp_int).max().getAsDouble());


            double[] emp_int2 = new double[iter];
            id= 0;
            intersection = new HashMap<>();
            map1 = frequencyMap_2bit[0].get(0);
            map2 = frequencyMap_2bit[1].get(0);

            for (Map.Entry<List<Double>, Integer> entry : map1.entrySet()) {
                List<Double> key = entry.getKey();
                if (map2.containsKey(key)) {
                    intersection.put(key, Math.min(entry.getValue(), map2.get(key)));
                }
            }
            for (Map.Entry<List<Double>, Integer> entry : intersection.entrySet()) {

                List<Double> listKey = entry.getKey();
                if(frequencyMap_2bit[0].get(0).get(listKey) < 1 || frequencyMap_2bit[1].get(0).get(listKey) < 1) {
                    continue;
                }
                double[] cp1 = clopperPearson((int)frequencyMap_2bit[0].get(0).get(listKey), iter, 0.05);
                double[] cp2 = clopperPearson((int)frequencyMap_2bit[1].get(0).get(listKey), iter, 0.05);
                double log_m1 = Math.abs(Math.log(cp1[0] / cp2[1]));
                double log_m2 = Math.abs(Math.log(cp1[1] / cp2[0]));
                emp_int2[id] = Math.abs(Math.max(log_m1, log_m2));
                if(Double.isInfinite(emp_int2[id])) {
                    emp_int[id] = 0;
                }
                id = id +1;
            }
            System.out.println("2 bit: " + Arrays.stream(emp_int2).max().getAsDouble());


        }
    }

    //限制数据取值范围也是[-1, 1] bin的个数和kv中的设置相同
    public double[] cal_boundary(int b1, int b2, int binNum){
        double step = (double)(b2 - b1) / binNum;
        double[] boundary = new double[2 * binNum];

        boundary[0] = -1;
        boundary[binNum*2-1] = 1;

        for(int i = 1; i < binNum*2 - 1; i = i + 2){
            boundary[i] = boundary[i - 1] + step;
            boundary[i + 1] = boundary[i];
        }
        return boundary;
    }

    public double laplace(double epsilon) {
        Random random = new Random();
        double sensitivity = 2.0;
        double scale = sensitivity / epsilon;
        // 生成拉普拉斯噪声
        double u = random.nextDouble() - 0.5;
        return -scale * Math.signum(u) * Math.log(1 - 2 * Math.abs(u));
    }

    // 对数据进行扰动
    public double[] perturbation(int value, double epsilon, double[] boundary) {
        double[] encode = new double[binNum];
        // 编码，将特定位置置为 1
        for (int i = 0; i < binNum; i++) {
            encode[i] = (i == value) ? 1.0 : 0.0;
        }

        // 添加拉普拉斯噪声
        for (int i = 0; i < binNum; i++) {
            encode[i] = encode[i] + laplace(epsilon);
        }

        return encode;
    }

    // 实现SHE函数
    public double[] SHE(int value, int b1, int b2, double epsilon, double theta) {
        double[] boundary = cal_boundary(b1, b2, binNum);
        double[] perturb_she = perturbation(value, epsilon, boundary);
        return perturb_she;
    }

    // 实现THE函数
    public double[] THE(int value, int b1, int b2, double epsilon, double theta) {
        double[] boundary = cal_boundary(b1, b2, binNum);
        double[] perturb =  perturbation(value, epsilon, boundary);
        double[] perturb_the = new double[binNum];
        for(int i = 0; i < binNum; i++){
            if(perturb[i] > theta){
                perturb_the[i] = 1;
            }else
                perturb_the[i] = 0;
        }
        return perturb_the;  //support
    }

    public static double[] clopperPearson(int k, int n, double alpha) {
        double lowerBound, upperBound;


        BetaDistribution betaLower = new BetaDistribution(k, n - k + 1);
        lowerBound = betaLower.inverseCumulativeProbability(alpha / 2);


        BetaDistribution betaUpper = new BetaDistribution(k + 1, n - k);
        upperBound = betaUpper.inverseCumulativeProbability(1 - alpha / 2);
        return new double[]{lowerBound, upperBound};
    }
    public static void main(String[] args)  throws Exception {
        new the();

    }

}

