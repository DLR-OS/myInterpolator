//
// This file is part of myInterpolator - image void filling
// Copyright (C) 2022  Dirk 'jtk' Frommholz, DLR OS-SEC
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

#ifndef COMMON_PROGRESS_H
#define COMMON_PROGRESS_H

#include "coolmath.h"
#include "StopWatch.h"


namespace common {


    /**
     * CProgress implements a progress indicator. It provides a
     * textural description of the current operation as well as
     * the progress percentage and elapsed/remaining time.
     *
     * The output routine can be adopted to a particular UI system
     * by subclassing.
     *
     * @author Dirk "jtk" Frommholz
     * @date January 27, 2010
     */
    class CProgress {
        private:


            /**
             * The text describing the current operation in progress.
             */
            std::string m_operation;


            /**
             * The progress of the current operation in percent (range 0.0 ... 
             * 1.0).
             */
            double m_percentage;


            /**
             * The previous progress percentage.
             */
            double m_oldPercentage;


            /**
             * The progress quantization in steps (i.e. display progress if the
             * percentage change is more than this value).
             */
            double m_percentageSteps;


            /**
             * The elapsed time since the initialization of the progress
             * indicator.
             */
            CStopWatch m_elaTime;

        
            /**
             * The remaining time in seconds. Interpolated from the elapsed
             * time and the percentage. If invalid, the value will be
             * negative.
             */        
            double m_remainingSecs;


            /**
             * True if the current operation has been aborted by the user.
             */
            bool m_abort;


        protected:


            /**
             * NeedToVisualizeProgress returns true if the progress needs to
             * be visualized because it has changed significantly.
             *
             */
            bool needToVisualizeProgress() {
                if (fp_equal(m_percentage, 0.0) || 
                fp_greater(m_percentage, 1.0) ||
                (fabs(m_oldPercentage-getPercentage())>m_percentageSteps)) {
                    m_oldPercentage=m_percentage;
                    return true;
                }

                return false;
            }
        

        public:


            /**
             * The constructor.
             */
            CProgress():m_operation(""), m_abort(false), 
            m_percentage(0.0), m_oldPercentage(0.0), m_remainingSecs(-1), 
            m_percentageSteps(0.05) { }


            /**
             * Init initializes the progress indicator using the given 
             * operation name. The percentage will be set to zero, and the 
             * timers will be reset too.
             *
             * @param f_operation the name of the operation to be started
             */
            void init(const std::string& f_operation) {
                m_operation=f_operation;
                m_percentage=0.0;
                m_oldPercentage=0.0;
                m_elaTime.start();
            }


            /**
             * GetElapsedTimeStr returns the elapsed time since the last 
             * call to init as a human readable string.
             *
             * @return the elapsed time of the current operation as a human
             * readable string
             */
            std::string getElapsedTimeStr() {
                m_elaTime.stop();
                return CStopWatch::secsToHr(m_elaTime.getPassedTime());
            }


            /**
             * GetRemainingTimeStr returns a human-readable string containing
             * the time the current operation will probably need to be
             * completely executed. The estimate is computed from the
             * time elapsed and the percentage.
             *
             */
            std::string getRemainingTimeStr() {
                m_elaTime.stop();            
                const double elaTime=m_elaTime.getPassedTime();            

                // cannot estimate very short remaining times
                if ((elaTime<1) || (getPercentage()<0.01)) {
                    return "n/a";
                }

                const double etaTime=elaTime/getPercentage();
                return CStopWatch::secsToHr(etaTime-elaTime);
            }


            /** 
             * SetOperation sets the descriptive text of the operation currently
             * in progress without changing the progress percentage of timers.
             *
             * @param f_operation the operation currently in progress
             */
            void setOperation(std::string f_operation) {
                m_operation=f_operation;
            }


            /**
             * GetOperation returns the operation currently in progress.
             *
             * @return the operation currently in progress
             */
            const std::string& getOperation() const {
                return m_operation;
            }        


            /**
             * SetPercentage sets the progress of the current operation in 
             * percent and recomputes the time values.
             *
             * @param f_percentage the percentage ranging from 0.0 ... 1.0
             */
            void setPercentage(double f_percentage) {
                m_percentage=f_percentage;
            }


            /**
             * SetPercentageSteps sets the percentage the progress must change
             * before it is visualized.
             *
             * @param f_percentageSteps the percentage quantization (i.e.
             * 0.05 -> display progress if it has changed more than 5%)
             */
            void setPercentageSteps(double f_percentageSteps) {
                m_percentageSteps=f_percentageSteps;
            }


            /**
             * GetPercentage returns the progress in percents (0.0 ... 1.0)
             * of the current operation.
             *
             * @return the progress percentage of the current operation
             */
            double getPercentage() const {
                return m_percentage;
            }


            /**
             * SetAbort sets the abort flag if the user wishes to abort the current
             * operation.
             *
             * @param f_abort true if the operation shall be aborted, false if not
             */
            void setAbort(bool f_abort) {
                m_abort=f_abort;
            }


            /**
             * GetAbort returns true if the user requests the current operation to be 
             * aborted.
             *
             * @return true on abort request, false otherwise
             */
            bool getAbort() {
                return m_abort;
            }

    
            /**
             * ShowProgress shows and updates the progress indicator using the 
             * underlying user interface. Therefore, this function has to be 
             * implemented by subclassing.
             *
             * @param f_force show the progress no matter if it's time to
             */
            virtual void showProgress(bool f_force=false) =0;


            /**
             * The destructor.
             */
            virtual ~CProgress() { }
    };

} // namespace common

#endif
