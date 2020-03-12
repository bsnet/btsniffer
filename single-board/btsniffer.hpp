/* 
 * Copyright 2019-2020 Marco Cominelli
 * Copyright 2019-2020 Francesco Gringoli
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BTDECODER_HPP
#define BTDECODER_HPP

/*
 * 8 CHANNELS
 */
#ifdef CHAN_8

const unsigned int srate = 1;
const unsigned int raw_srate = 8;

#define FILTER_TAP_NUM 17
static double filter_taps[FILTER_TAP_NUM] = {
  0.006024041056991322,
  0.018694160465501236,
  0.011702650731212858,
  -0.023756813610156347,
  -0.0583711817341587,
  -0.022338894580199133,
  0.1140137071558402,
  0.2831550731366491,
  0.36044780676612975,
  0.2831550731366491,
  0.1140137071558402,
  -0.022338894580199133,
  -0.0583711817341587,
  -0.023756813610156347,
  0.011702650731212858,
  0.018694160465501236,
  0.006024041056991322
};

#endif // CHAN_8


/*
 * 16 CHANNELS
 */
#ifdef CHAN_16

const unsigned int srate = 1;
const unsigned int raw_srate = 16;

#define FILTER_TAP_NUM 33
static double filter_taps[FILTER_TAP_NUM] =
{
  -0.007876653238814611,
  -0.0045078023642330135,
  -0.004647176836099708,
  -0.0036762062668577035,
  -0.001232789068763745,
  0.002974057039622697,
  0.009113836327699427,
  0.017189231710653903,
  0.027012802038788086,
  0.03820205049740163,
  0.0501884329045113,
  0.06225825641438132,
  0.07362529190934415,
  0.08349929494101796,
  0.09114736692124643,
  0.09597796817945085,
  0.09764686418550765,
  0.09597796817945085,
  0.09114736692124643,
  0.08349929494101796,
  0.07362529190934415,
  0.06225825641438132,
  0.0501884329045113,
  0.03820205049740163,
  0.027012802038788086,
  0.017189231710653903,
  0.009113836327699427,
  0.002974057039622697,
  -0.001232789068763745,
  -0.0036762062668577035,
  -0.004647176836099708,
  -0.0045078023642330135,
  -0.007876653238814611
};

#endif // CHAN_16


#endif // BTDECODER_HPP
